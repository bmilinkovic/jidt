package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import infodynamics.utils.MatrixUtils;
import junit.framework.TestCase;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Pedro A.M. Mediano (<a href="pmediano at pm.me">email</a>,
 * <a href="http://www.doc.ic.ac.uk/~pam213">www</a>)
 */
public class EffectiveInformationDiscreteTester extends TestCase {

    protected RandomGenerator rg = new RandomGenerator();

    public void testAddObservations() {

        int[] t0 = { 0, 0, 0, 0};
        int[] t1 = { 1, 1, 1, 0};
        int[] t2 = { 1, 1, 1, 0};
        int[] t3 = { 0, 0, 1, 0};
        int[] t4 = { 1, 1, 1, 0};
        int[] t5 = { 0, 0, 1, 0};
        int[][] states0 = {t0, t1, t2, t3, t4, t5};
        int base = 2, dims = 4;
        boolean exceptionCaught = false;
        EffectiveInformationCalculatorDiscrete eicd;
        eicd = new EffectiveInformationCalculatorDiscrete(base, dims);
        try {
          eicd.initialise();
          eicd.setObservations(states0);
        } catch (Throwable e) {
          exceptionCaught = true;
        }
        assertFalse(exceptionCaught);

        int[][] states1 = eicd.getData();
        assertTrue(Arrays.equals(states0, states1));
    }

    public void testGenerativeModel() {
        // Example adapted from Wikipedia article on IIT.
        // https://en.wikipedia.org/wiki/Integrated_information_theory

        // Generate model.
        RandomGenerator rg = new RandomGenerator();
        int duration = 1000000;

        int[] var0a = rg.generateRandomInts(duration, 2);
        var0a[0] = 0;
        int[] var1a = new int[duration];
        var1a[0] = 0;
        System.arraycopy(var0a, 0, var1a, 1, duration - 1);
        int[][] states = {var0a, var1a};
        states = MatrixUtils.transpose(states);

        boolean exceptionCaught = false;
        int base = 2, dims = 2;
        EffectiveInformationCalculatorDiscrete eicd;
        eicd = new EffectiveInformationCalculatorDiscrete(base, dims);
        try {
          eicd.initialise();
          eicd.setObservations(states);
        } catch (Throwable e) {
          exceptionCaught = true;
        }
        assertFalse(exceptionCaught);

        // The MI and effective information for this model should be very close
        // to 1, but there is a margin of error that needs to be considered.
        double marginOfError = 0.001;

        // Calculate mutual information for system.
        double system = eicd.computeForSystem();
        double diffMi = Math.abs(system - 1);
        assertTrue(diffMi < marginOfError);

        // Calculate EI for the partition.
        List<List<Integer>> p = new ArrayList<List<Integer>>();
        p.add(new ArrayList<>(Arrays.asList(0)));
        p.add(new ArrayList<>(Arrays.asList(1)));
        double ei = eicd.computeForPartition(p);
        double diffEi = Math.abs(ei - 1);
        assertTrue(diffEi < marginOfError);

    }

    public void testUnnormalisedPDF() {
        // Test that calculator complains when given an unnormalised PDF
        double[][] pdf = new double[][] {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
        EffectiveInformationCalculatorDiscrete eicd;
        eicd = new EffectiveInformationCalculatorDiscrete(2, 2);
        boolean exceptionCaught = false;
        try {
          eicd.setJointPDF(pdf);
        } catch (Throwable e) {
          exceptionCaught = true;
        }
        assertTrue(exceptionCaught);
    }

    public void testComputePartition() throws Exception {

        double tol = 0.001;

        EffectiveInformationCalculatorDiscrete eicd = new EffectiveInformationCalculatorDiscrete(2, 4);
        List<List<Integer>> pl = new ArrayList<>(List.of(new ArrayList<>(List.of(0)), new ArrayList<>(List.of(1,2,3))));
        int[] pv = new int[] {0, 1, 1, 1};

        int[][] data = rg.generateRandomInts(100, 4, 2);
        eicd.initialise();
        eicd.setObservations(data);
        assertEquals(eicd.computeForPartition(pl), eicd.computeForPartition(pv), tol);
    }

}
