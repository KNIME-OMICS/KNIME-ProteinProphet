package de.mpc.tools.knimeproteinprophet;

import org.knime.base.node.util.exttool.ExtToolStderrNodeView;
import org.knime.base.node.util.exttool.ExtToolStdoutNodeView;
import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ProteinProphet" Node.
 * KNIME node to perform ProteinProphet inference
 *
 * @author julianu
 */
public class ProteinProphetNodeFactory 
        extends NodeFactory<ProteinProphetNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ProteinProphetNodeModel createNodeModel() {
        return new ProteinProphetNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 2;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<ProteinProphetNodeModel> createNodeView(final int viewIndex,
            final ProteinProphetNodeModel nodeModel) {
    	if (viewIndex == 0) {
    		return new ExtToolStdoutNodeView<ProteinProphetNodeModel>(nodeModel);
    	} else if (viewIndex == 1) {
    		return new ExtToolStderrNodeView<ProteinProphetNodeModel>(nodeModel);
    	}
    	return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new ProteinProphetNodeDialog();
    }

}

