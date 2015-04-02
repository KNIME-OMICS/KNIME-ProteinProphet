package de.mpc.tools.knimeproteinprophet;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;


public class ProteinProphetRunnable implements Runnable {
	
	/** the input pepXML files */
	private List<String> pepXMLfiles;
	
	/** the FASTA database */
	private String fastaFile;
	
	/** the used enzyme (in pepXML annotation, i.e. T=Trypsin,  C=Chymotrypsin etc.)*/
	private String enzyme;
	
	/** the minimal peptide probability to use */
	private Double peptide_prob;
	
	/** the decoy prefix */
	private String decoyPrefix;
	
	/** maximal number of threads to use */
	private Integer threads;
	
	/** path to the executable of xinteract */
	private String execXinteract;
	
	/** path to the executable of ProteinProphet */
	private String execProteinProphet;
	
	/** path to the (temporary) directory */
	private String executionDirectory;
	
	/** the STD output*/
	private List<String> output;
	
	/** the STDERR output */
	private List<String> errorOutput;
	
	/** the currently running thread */
	private Process runningProcess = null;
	
	/** the final protXML file */
	private String protXMLFile;
	
	/** the final tab separated file */
	private String excelFile;
	
	
	public ProteinProphetRunnable(List<String> pepXMLfiles, String fastaFile, String enzyme, Double peptide_prob, String decoyPrefix,
			Integer threads, String execXinteract, String execProteinProphet, String executionDirectory,
			List<String> output, List<String> errorOutput) {
		this.pepXMLfiles = pepXMLfiles;
		this.fastaFile = fastaFile;
		this.enzyme = enzyme;
		this.peptide_prob = peptide_prob;
		this.decoyPrefix = decoyPrefix;
		this.threads = threads;
		this.execXinteract = execXinteract;
		this.execProteinProphet = execProteinProphet;
		this.executionDirectory = executionDirectory;
		this.output = output;
		this.errorOutput = errorOutput;
		this.protXMLFile = null;
		this.excelFile = null;
    }
    
    
	@Override
	public void run() {
		ProcessBuilder processB = new ProcessBuilder(
				execXinteract,
				"-D" + fastaFile,
				"-e" + enzyme,
				"-nP",
				"-Ot",
				"-d" + decoyPrefix,
				"-THREADS=" + threads,
				"-i",
				"-N" + executionDirectory + File.separator + "xinteractout.pep.xml");
		processB.command().addAll(pepXMLfiles);
		
		try {
			BufferedReader stdOut = null;
			BufferedReader stdError = null;
			
			runningProcess = processB.start();
			
			// read the output from xinteract
			if (runningProcess != null) {
				stdOut = new BufferedReader(new InputStreamReader(runningProcess.getInputStream()));
				stdError = new BufferedReader(new InputStreamReader(runningProcess.getErrorStream()));
			}
			
			String s;
			while ((runningProcess != null) && ((s = stdOut.readLine()) != null)) {
				output.add(s);
			}
			while ((runningProcess != null) && ((s = stdError.readLine()) != null)) {
				errorOutput.add(s);
			}
			
			runningProcess = null;
			stdOut = null;
			stdError = null;
			processB = new ProcessBuilder(
					execProteinProphet,
					executionDirectory + File.separator + "xinteractout.ipro.pep.xml",
					executionDirectory + File.separator + "proteinprophet.protXML",
					"IPROPHET",
					"MINPROB" + peptide_prob,
					"NOPLOT",
					"EXCELPEPS");
			
			runningProcess = processB.start();
			
			if (runningProcess != null) {
				stdOut = new BufferedReader(new InputStreamReader(runningProcess.getInputStream()));
				stdError = new BufferedReader(new InputStreamReader(runningProcess.getErrorStream()));
			}
			
			// read the output from ProteinProphet
			while ((runningProcess != null) && ((s = stdOut.readLine()) != null)) {
				output.add(s);
			}
			while ((runningProcess != null) && ((s = stdError.readLine()) != null)) {
				errorOutput.add(s);
			}
			
			protXMLFile = executionDirectory + File.separator + "proteinprophet.protXML";
			excelFile = executionDirectory + File.separator + "proteinprophet.xls";
		} catch (IOException e) {
			ProteinProphetNodeModel.logger.error("Error while executing", e);
			protXMLFile = null;
			excelFile = null;
		}
	}
	
	
	/**
	 * returns the path to the created protXML file
	 * @return
	 */
	public String getProtXMLFile() {
		return protXMLFile;
	}
	
	
	/**
	 * returns the path to the created excel file
	 * @return
	 */
	public String getExcelFile() {
		return excelFile;
	}
}
