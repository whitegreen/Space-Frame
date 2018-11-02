package test;

import java.text.DecimalFormat;

import frame.Frame2D;

// A first course in the Finite Element Method 5th ed, Logan, chapter.5 Frame
//Hao Hua, Southeast University, whitegreen@163.com

public class TestFrame2D {
	DecimalFormat df = new DecimalFormat("###.##");

	public static void main(String[] args) {
		new TestFrame2D();
	}

	public TestFrame2D() {
		//Logan p240
//		double[][] nodes = new double[][] { { 0, 0 }, { 0, 120 }, { 120, 120 }, { 120, 0 } }; // inch
//		int[][] members = new int[][] { { 0, 1 }, { 1, 2 }, { 2, 3 } }; // node ids
//		double[][] loads = new double[nodes.length][];
//		loads[1] = new double[] { 10000, 0, 0 }; // lb
//		loads[2] = new double[] { 0, 0, 5000 };
//		boolean[][] constrained = new boolean[nodes.length][];
//		constrained[0] = new boolean[] { true, true, true }; // RX, RY, moment;
//		constrained[3] = new boolean[] { true, true, true };
//		double[] As = new double[members.length];
//		double[] Es = new double[members.length];
//		double[] Is = new double[] { 200, 100, 200 }; // in
//		for (int i = 0; i < members.length; i++) {
//			As[i] = 10; // in
//			Es[i] = 3E7; // psi
//		}

		// 246C distributed loads,
		double[][] nodes = new double[][] { { 0, 0 }, { 360, 360 }, { 840, 360 } }; // in
		int[][] members = new int[][] { { 0, 1 }, { 1, 2 } }; // node ids
		Double[] memb_selfweis = new Double[members.length]; //self-weight of beam, unit: N or lb.
		memb_selfweis[0] = 6363.96; // self-weight of beam, unit: N or lb.
		memb_selfweis[1] = 2400.00; 
		double[][] loads = new double[nodes.length][3];
		boolean[][] constrained = new boolean[nodes.length][];
		constrained[0] = new boolean[] { true, true, true }; // RX, RY, moment;
		constrained[2] = new boolean[] { true, true, true };
		double[] As = new double[] { 100, 140 };
		double[] Es = new double[] { 30000, 35000 };
		double[] Is = new double[] { 1000, 800 }; //principal moment of inertia

		//Frame2D fr = new Frame2D(nodes, members, loads, constrained, As, Es, Is, null); //no self weight
		Frame2D fr = new Frame2D(nodes, members, loads, constrained, As, Es, Is, memb_selfweis);
		
		println("nodal displace / rotation: ");
		for (int i = 0; i < fr.lenA; i++) 
			println(fr.symbols[i] + "   " + df.format(fr.UA[i]));

		System.out.println("member (local) forces / bending moments");
		for (int i = 0; i < fr.member_forces.length; i++) {
			double[] a = fr.member_forces[i];
			System.out.print("member " + i + ": ");
			for (int j = 0; j < a.length; j++) {
				if (Frame2D.DIM == j) {
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
