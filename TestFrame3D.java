package test;

import java.text.DecimalFormat;
import frame.Frame3D;

/**
 * A first course in the Finite Element Method 5th ed, Logan, chapter.5 Frame
 * 
 * @author Hao Hua, Southeast University, whitegreen@163.com
 *
 **/

public class TestFrame3D {
	private DecimalFormat df = new DecimalFormat("###.##");

	public static void main(String[] args) {
		new TestFrame3D();
	}

	public TestFrame3D() {
		double[][] nodes = new double[][] { { 53, 87, 34 }, { 59, 71, 60 }, { 42, 35, 44 }, { 22, 0, 13 }, { 0, 14, 66 }, { -5, -2, 0 } };
		int[][] members = new int[][] { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 2, 4 }, { 3, 5 } }; // node ids
		double[][] loads = new double[nodes.length][6];   // nodal: fx, fy, fz, mx, my, mz
		loads[1] = new double[] { 100, -55, 99, 1000, -2000, 900 }; // fx, fy, fz, mx, my, mz
		loads[3] = new double[] { -50, 100, 0,  800, 1200, 1000 };
		Double[] memb_selfweis = new Double[members.length];  //self-weight of beam, unit: N or lb.
		memb_selfweis[0] = 311.1269837;
		memb_selfweis[2] = 508.5272854;  //always positive
		boolean[][] constrained = new boolean[nodes.length][];
		constrained[0] = new boolean[] { true, true, true, false, true, false }; // RX, RY, RZ, mx, my, mz;
		constrained[4] = new boolean[] { true, true, false, true, true, true };
		constrained[5] = new boolean[] { true, false, true, true, false, false };
		double[] As = new double[] { 10, 8, 6, 12, 14 };
		double[] Es = new double[] { 30000, 40000, 50000, 10000, 20000, };
		double[] Gs = new double[] { 8000, 10000, 12000, 16000, 14000 }; // shear modulus, or modulus of rigidity
		double[] Js = new double[] { 100, 80, 60, 30, 70 }; // torsion constant
		double[] Iy = new double[] { 150, 100, 60, 120, 90 };  // principal moment of inertia
		double[] Iz = new double[] { 100, 80, 70, 130, 120 };  // principal moment of inertia
		
		Frame3D fr = new Frame3D(nodes, members, loads, constrained);
		fr.analyzeFrame(As, Es, Gs, Js, Iy, Iz, memb_selfweis);

		println("nodal displace / rotation: ");
		for (int i = 0; i < fr.lenA; i++) 
			println(fr.symbols[i] + "   " + df.format( fr.UA[i]));

		println("member (local) forces / bending moments");
		for (int i = 0; i < fr.member_forces.length; i++) {
			double[] a = fr.member_forces[i];
			System.out.print("member " + i + ": ");
			for (int j = 0; j < a.length; j++) {
				if (6 == j) {
					println("");
					System.out.print("         ");
				}
				System.out.print(df.format(a[j]) + ",");
			}
			println("");
		}
	}
	private void println(String s){
		System.out.println(s);
	}

}
