package de.mpc.tools.knimeproteinprophet;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.knime.base.node.util.exttool.ExtToolOutputNodeModel;
import org.knime.core.data.uri.IURIPortObject;
import org.knime.core.data.uri.URIContent;
import org.knime.core.data.uri.URIPortObject;
import org.knime.core.data.uri.URIPortObjectSpec;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortObjectSpec;
import org.knime.core.node.port.PortType;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;


/**
 * This is the model implementation of ProteinProphet.
 * KNIME node to perform ProteinProphet inference
 *
 * @author julianu
 */
public class ProteinProphetNodeModel extends ExtToolOutputNodeModel {
    
	// the logger instance
	protected static final NodeLogger logger = NodeLogger
	        .getLogger(ProteinProphetNodeModel.class);
	    
    
	static final Map<String, String> enzymeNameToShort;
	static {
		enzymeNameToShort = new HashMap<String, String>();
		
		enzymeNameToShort.put("Trypsin", "T");
		enzymeNameToShort.put("StrictTrypsin", "S");
		enzymeNameToShort.put("Chymotrypsin", "C");
		enzymeNameToShort.put("RalphTrypsin", "R");
		enzymeNameToShort.put("AspN", "A");
		enzymeNameToShort.put("GluC", "G");
		enzymeNameToShort.put("GluC Bicarb", "B");
		enzymeNameToShort.put("CNBr", "M");
		enzymeNameToShort.put("Trypsin/CNBr", "D");
		enzymeNameToShort.put("Chymotrypsin/AspN/Trypsin", "3");
		enzymeNameToShort.put("Elastase", "E");
		enzymeNameToShort.put("LysC / Trypsin_K (cuts after K not before P)]", "K");
		enzymeNameToShort.put("LysN (cuts before K)]", "L");
		enzymeNameToShort.put("LysN Promisc (cuts before KASR)]", "P");
		enzymeNameToShort.put("Nonspecific or None", "N");
	}
	
	static final String CFGKEY_ENZYME = "Enzyme";
	static final String[] ALLOWED_ENZYMES = enzymeNameToShort.keySet().toArray(new String[1]);
	static final String DEFAULT_ENZYME = ALLOWED_ENZYMES[0];
	
	
	private final SettingsModelString m_enzyme =
			new SettingsModelString(ProteinProphetNodeModel.CFGKEY_ENZYME, ProteinProphetNodeModel.DEFAULT_ENZYME);
	
	
	static final String CFGKEY_DECOYPREFIX = "Decoyprefix";
	static final String DEFAULT_DECOYPREFIX = "decoy_";
	
	private final SettingsModelString m_decoyprefix =
			new SettingsModelString(ProteinProphetNodeModel.CFGKEY_DECOYPREFIX, ProteinProphetNodeModel.DEFAULT_DECOYPREFIX);
	
	
	static final String CFGKEY_THREADS = "Threads";
	static final Integer DEFAULT_THREADS = 1;
	
	private final SettingsModelInteger m_threads =
			new SettingsModelInteger(ProteinProphetNodeModel.CFGKEY_THREADS, ProteinProphetNodeModel.DEFAULT_THREADS);
	
	
	/** the executable for xinteract */
	private File execXinteract = null;
	
	/** the executable for ProteinProphet */
	private File execProteinProphet = null;
	
	
	/** the actual execution thread */
	private Thread executionThread = null;
	
	
	/**
	 * Constructor for the node model.
	 */
	protected ProteinProphetNodeModel() {
		// two incoming URI ports, one outgoing URI port
		super(new PortType[]{IURIPortObject.TYPE, IURIPortObject.TYPE},
				new PortType[]{IURIPortObject.TYPE, IURIPortObject.TYPE});
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObject[] execute(PortObject[] inObjects, ExecutionContext execContext)
			throws Exception {
		
		// create a working directory
		Random randomNumberGenerator = new Random();
		int num = randomNumberGenerator.nextInt(Integer.MAX_VALUE);
		File dir = new File(System.getProperty("java.io.tmpdir") + File.separator
				+ String.format("%s%06d", "PPinference", num));
		
		while (dir.exists()) {
			num = randomNumberGenerator.nextInt(Integer.MAX_VALUE);
			dir = new File(System.getProperty("java.io.tmpdir") + File.separator
					+ String.format("%s%06d", "PPinference", num));
		}
		dir.mkdirs();
		dir.deleteOnExit();
		
		// get the input pepXML files
		IURIPortObject filesPort = (IURIPortObject) inObjects[0];
		List<URIContent> uris = filesPort.getURIContents();
		List<String> inputFiles = new ArrayList<String>(uris.size()); 
		for (URIContent uric : uris) {
            URI uri = uric.getURI();
            inputFiles.add(new File(uri).getAbsolutePath());
        }
		
		// get the input FASTA files
		filesPort = (IURIPortObject) inObjects[1];
		uris = filesPort.getURIContents();
		String fastaFile = null; 
		for (URIContent uric : uris) {
            URI uri = uric.getURI();
            fastaFile = new File(uri).getAbsolutePath();
            break;
        }
		
		
		// create correct enzyme
		String enzyme = enzymeNameToShort.get(m_enzyme.getStringValue());
		
		// check the input files, whether it has the enzyme in the "<msms_run_summary>" tag
		ListIterator<String> inputFilesIt = inputFiles.listIterator();
		while (inputFilesIt.hasNext()) {
			String file = inputFilesIt.next();
			String newFile = checkInputFile(file, enzyme, dir);
			
			if (!newFile.equals(file)) {
				inputFilesIt.set(newFile);
			}
		}
		
		
		LinkedList<String> externalOutput = new LinkedList<String>();
		LinkedList<String> externalErrorOutput = new LinkedList<String>();
		
		ProteinProphetRunnable pprunner =
				new ProteinProphetRunnable(inputFiles, fastaFile, enzyme, m_decoyprefix.getStringValue(), m_threads.getIntValue(),
						execXinteract.getAbsolutePath(), execProteinProphet.getAbsolutePath(), dir.getAbsolutePath(),
						externalOutput, externalErrorOutput);
		
		executionThread = new Thread(pprunner);
		executionThread.start();
		
		while (executionThread.isAlive()) {
			executionThread.join(1000);
			
			try {
				execContext.checkCanceled();
			} catch (CanceledExecutionException ex) {
				executionThread.interrupt();
				executionThread.join();
			}
		}
		
		
		String protXMLFile = pprunner.getProtXMLFile();
		List<URIContent> outProtXML = new ArrayList<URIContent>();
		List<URIContent> outXLS = new ArrayList<URIContent>();
		
		if ((protXMLFile != null) && Files.exists(new File(protXMLFile).toPath(), new LinkOption[]{})) {
			outProtXML.add(new URIContent(new File(protXMLFile).toURI(), "protXML"));
			outXLS.add(new URIContent(new File(pprunner.getExcelFile()).toURI(), "xls"));
			
			setExternalOutput(externalOutput);
			setExternalErrorOutput(externalErrorOutput);
		} else {
			setFailedExternalOutput(externalOutput);
			setFailedExternalErrorOutput(externalErrorOutput);
			throw new Exception("Error while executing ProteinProphet.");
		}
		
		URIPortObject outProtXMLPort = new URIPortObject(outProtXML);
		URIPortObject outXLSPort = new URIPortObject(outXLS);
		
        return new PortObject[]{outProtXMLPort, outXLSPort};
	}
	
	
	/**
	 * Checks the input file for errors and corrects them while copying the
	 * file and returning the new filename.<p>
	 * Additionally, it needs to change the search-engine name to avoid trouble
	 * with unneeded corrections.
	 * 
	 * @param fileName
	 * @return
	 */
	private String checkInputFile(String fileName, String enzymeShort, File tmpDir) throws FileNotFoundException, IOException{
		String newFilename = fileName;
		
		boolean containsEnzyme = false;
		
		// pre-check the file
		InputStream fis = new FileInputStream(fileName);
		InputStreamReader isr = new InputStreamReader(fis);
		BufferedReader br = new BufferedReader(isr);
		
		String line;
		while ((line = br.readLine()) != null) {
			
			if (line.contains("<sample_enzyme")) {
				containsEnzyme = true;
			}
		}
		br.close();
		
		
		// now copy the file and perform corrections
		fis = new FileInputStream(fileName);
		isr = new InputStreamReader(fis);
		br = new BufferedReader(isr);
		
		newFilename = tmpDir.getAbsolutePath() + File.separator + new File(fileName).getName();
		FileWriter out = new FileWriter(newFilename);
		BufferedWriter bw = new BufferedWriter(out);
		
		Pattern se_pattern = Pattern.compile(".*search_engine=\"([^\"]*)\".*");
		
		while ((line = br.readLine()) != null) {
			Matcher se_matcher = se_pattern.matcher(line);
			if (se_matcher.matches()) {
				line = line.replaceAll("search_engine=\"[^\"]*\"", "search_engine=\"" + se_matcher.group(1) + "-correct\"");
			}
			
			bw.append(line);
			bw.append("\n");
			
			if (!containsEnzyme && line.contains("<msms_run_summary")) {
				// add the enzyme tag here
				logger.warn(fileName + " needs to add the enzyme tag.");
				bw.append(createEnzymeTag(enzymeShort));
				bw.append("\n");
			}
		}
		
		bw.close();
		br.close();
		
		return newFilename;
	}
	
	
	/**
	 * create the tag for the enzyme in pepXML
	 * @param enzyme
	 * @return
	 */
	private String createEnzymeTag(String enzyme) {
		String name = null;
		String cut = null;
		String no_cut = null;
		String sense = null;
		
		if (enzyme == "T") {
			name = "trypsin";
			cut = "KR";
			no_cut = "P";
			sense = "C";
		}
		
		return "\t<sample_enzyme name=\"" + name + "\">\n\t\t" +
				"<specificity cut=\"" + cut + "\" no_cut=\"" + no_cut + "\" sense=\"" + sense + "\"/>\n\t" +
				"</sample_enzyme>";
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected PortObjectSpec[] configure(final PortObjectSpec[] inSpecs)
			throws InvalidSettingsException {
		
		// check for the executables
		String path;
		try {
			path = this.getClass().getProtectionDomain().getCodeSource().getLocation().toURI().getPath();
			if (!path.endsWith(File.separator)) {
				// we are in the jar, only get the path to it
				path = path.substring(0, path.lastIndexOf(File.separator) + 1);
			}
			
			path += "executables" + File.separator;
			
			execXinteract = new File(path + File.separator + "xinteract");
			if (!Files.exists(execXinteract.toPath(), new LinkOption[]{})) {
				throw new InvalidSettingsException("Failed to find matching binary for xinteract in '" + path + "'.");
			}
			
			execProteinProphet = new File(path + File.separator + "ProteinProphet");
			if (!Files.exists(execProteinProphet.toPath(), new LinkOption[]{})) {
				throw new InvalidSettingsException("Failed to find matching binary for ProteinProphet in '" + path + "'.");
			}
		} catch (URISyntaxException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
        PortObjectSpec[] out_spec = new PortObjectSpec[2];
        out_spec[0] = new URIPortObjectSpec(new String[]{"protxml", "protXML", "XML"});
        out_spec[1] = new URIPortObjectSpec(new String[]{"xls"});
        
        return out_spec;
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		m_enzyme.saveSettingsTo(settings);
		m_decoyprefix.saveSettingsTo(settings);
		m_threads.saveSettingsTo(settings);
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_enzyme.loadSettingsFrom(settings);
		m_decoyprefix.loadSettingsFrom(settings);
		m_threads.loadSettingsFrom(settings);
	}
	
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_enzyme.validateSettings(settings);
		m_decoyprefix.validateSettings(settings);
		m_threads.validateSettings(settings);
	}
}
