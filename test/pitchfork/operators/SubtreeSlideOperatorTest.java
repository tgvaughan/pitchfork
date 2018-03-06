package pitchfork.operators;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.Assert;
import org.junit.Test;

public class SubtreeSlideOperatorTest {

    /* Tests WITHOUT coalescence node attachment */

    @Test
    public void testComputeYoungerAttachmentPointProb1() {

        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 0.0,
                "weight", 1.0);

        SubtreeSlideOperator.AttachmentPoint ap = stsOp.new AttachmentPoint();
        Node edgeParentNode = tree.getRoot();
        ap.attachmentEdgeBase = tree.getNode(0);
        ap.attachmentHeight = 1.2;

        stsOp.computeYoungerAttachmentPointProb(ap, edgeParentNode, stsOp.getCurrentLambda());

        Assert.assertEquals(-3.894639484342174, ap.logProb, 1e-10);
    }

    @Test
    public void testComputeOlderAttchmentPointProb1() {
        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 0.0,
                "weight", 1.0);

        SubtreeSlideOperator.AttachmentPoint ap = stsOp.new AttachmentPoint();
        Node edgeParentNode = tree.getNode(0).getParent();
        ap.attachmentEdgeBase = tree.getRoot();
        ap.attachmentHeight = 3.7;

        stsOp.computeOlderAttachmentPointProb(ap, edgeParentNode, stsOp.getCurrentLambda());

        Assert.assertEquals(-2.979270081560007, ap.logProb, 1e-10);
    }

    @Test
    public void testGetYoungerAttachmentPoint1() {
        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 0.0,
                "weight", 1.0);

        Node edgeBaseNode = tree.getNode(2);
        Node edgeBaseParent = edgeBaseNode.getParent();
        double lambda = stsOp.getCurrentLambda();

        SubtreeSlideOperator.AttachmentPoint ap = null;

        boolean succeeded = false;
        do {
            try {
                ap = stsOp.getYoungerAttachmentPoint(edgeBaseNode, edgeBaseParent, lambda);
                succeeded = true;
            } catch (SubtreeSlideOperator.AttachmentException ex) {
            }
        } while(!succeeded);

        System.out.println(ap);

        double logP1 = ap.logProb;

        ap.logProb = 0.0;
        stsOp.computeYoungerAttachmentPointProb(ap, edgeBaseParent, lambda);

        double logP2 = ap.logProb;

        Assert.assertEquals(logP1, logP2, 1e-10);
    }

    @Test
    public void testGetOlderAttachmentPoint1() {
        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 0.0,
                "weight", 1.0);

        Node startNode = tree.getNode(0).getParent();
        double lambda = stsOp.getCurrentLambda();

        SubtreeSlideOperator.AttachmentPoint ap =
                stsOp.getOlderAttachmentPoint(startNode, stsOp.getCurrentLambda());

        System.out.println(ap);

        double logP1 = ap.logProb;

        ap.logProb = 0.0;
        stsOp.computeOlderAttachmentPointProb(ap, startNode, lambda);

        double logP2 = ap.logProb;

        Assert.assertEquals(logP1, logP2, 1e-10);
    }

    /* Tests WITH coalescence node attachment */

    @Test
    public void testComputeOlderAttachmentPointProb2() {

        Tree tree = new TreeParser("((A:1,B:1):1,C:2):0.0");

        double pa = 0.5;

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", pa,
                "weight", 1.0);

        SubtreeSlideOperator.AttachmentPoint ap = stsOp.new AttachmentPoint();
        Node edgeBaseNode = tree.getNode(1);
        Node edgeParentNode = edgeBaseNode.getParent();
        ap.attachmentEdgeBase = tree.getRoot();
        ap.attachmentHeight = ap.attachmentEdgeBase.getHeight();

        stsOp.computeOlderAttachmentPointProb(ap, edgeParentNode, stsOp.getCurrentLambda());

        System.out.println(ap);

        Assert.assertEquals(Math.log(pa), ap.logProb, 1e-10);
    }

    @Test
    public void testGetYoungerAttachmentPoint2() {
        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 1.0,
                "weight", 1.0);

        Node edgeBaseNode = tree.getNode(2);
        Node edgeBaseParent = edgeBaseNode.getParent();
        double lambda = stsOp.getCurrentLambda();

        SubtreeSlideOperator.AttachmentPoint ap = null;

        boolean succeeded = false;
        do {
            try {
                ap = stsOp.getYoungerAttachmentPoint(edgeBaseNode, edgeBaseParent, lambda);
                succeeded = true;
            } catch (SubtreeSlideOperator.AttachmentException ex) {
            }
        } while(!succeeded);

        System.out.println(tree);
        System.out.println(ap);

        double logP1 = ap.logProb;

        ap.logProb = 0.0;
        stsOp.computeYoungerAttachmentPointProb(ap, edgeBaseParent, lambda);

        double logP2 = ap.logProb;

        Assert.assertEquals(logP1, logP2, 1e-10);
    }

    @Test
    public void testGetOlderAttachmentPoint2() {
        Tree tree = new TreeParser("((A:1,B:1):1,C:3):0.0");

        SubtreeSlideOperator stsOp = new SubtreeSlideOperator();
        stsOp.initByName(
                "tree",tree,
                "relSize",0.15,
                "probCoalAttach", 1.0,
                "weight", 1.0);

        Node startNode = tree.getNode(0).getParent();
        double lambda = stsOp.getCurrentLambda();

        SubtreeSlideOperator.AttachmentPoint ap =
                stsOp.getOlderAttachmentPoint(startNode, stsOp.getCurrentLambda());

        System.out.println(tree);
        System.out.println(ap);

        double logP1 = ap.logProb;

        ap.logProb = 0.0;
        stsOp.computeOlderAttachmentPointProb(ap, startNode, lambda);

        double logP2 = ap.logProb;

        Assert.assertEquals(logP1, logP2, 1e-10);
    }

}
