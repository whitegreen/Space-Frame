package test;

import java.text.DecimalFormat;
import java.util.ArrayList;

import peasy.PeasyCam;
import processing.core.PApplet;
import frame.Frame3D;
import frame.M;

public class Building2b extends PApplet {  // unit:m, a simple timber frame, a point load is added (option)
	private DecimalFormat df = new DecimalFormat("###.##");
	private static final double sec_w=0.15; //m
	private static final double wood_density=450;  // kg/m3
	private static final double earth_gravity= 9.8;// N/kg
	private double[][] nodes;
	private int[][] members;
	private boolean[][] constrained;
	private Frame3D fr;
	private double val_max, val_min;
	private double sc=50; //display
	private double[][] res;
	

	public void setup() {
		size(1000, 800, P3D);
		new PeasyCam(this, 200);

		double wid = 3.333; 
		double hei = 2.400;
		double high = 2.70;
		nodes = new double[6 * 3][];
		for (int k = 0; k < 3; k++)
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 3; j++)
					nodes[k * 6 + i * 3 + j] = new double[] { i * wid, (j - 1) * hei, -k * high };
		ArrayList<int[]> elist = new ArrayList<int[]>();
		for (int k = 0; k < 2; k++) {
			elist.add(new int[] { 6 * k + 0, 6 * k + 1 });
			elist.add(new int[] { 6 * k + 1, 6 * k + 2 });
			elist.add(new int[] { 6 * k + 3, 6 * k + 4 });
			elist.add(new int[] { 6 * k + 4, 6 * k + 5 });
			elist.add(new int[] { 6 * k + 0, 6 * k + 3 });
			elist.add(new int[] { 6 * k + 1, 6 * k + 4 });
			elist.add(new int[] { 6 * k + 2, 6 * k + 5 });
		}
		for (int i = 0; i < 6; i++) {
			elist.add(new int[] { i + 0, i + 6 });
			elist.add(new int[] { i + 6, i + 12 });
		}
		members = new int[elist.size()][];
		for (int i = 0; i < members.length; i++)
			members[i] = elist.get(i);
				
		double[][] loads = new double[nodes.length][6]; // nodal: fx, fy, fz, mx, my, mz
//		loads[1] = new double[] { 100, -55, 99, 1000, -2000, 900 }; // fx, fy, fz, mx, my, mz
//		loads[3] = new double[] { -50, 100, 0, 800, 1200, 1000 };
		Double[] memb_selfweis = new Double[members.length]; // self-weight of beam, unit: N
		for (int i = 0; i < members.length; i++) {
			int[] m = members[i];
			double[] pa = nodes[m[0]];
			double[] pb = nodes[m[1]];
			double len=M.dist(pa, pb);
			double vol= sec_w*sec_w * len;
			memb_selfweis[i]= earth_gravity*(wood_density*vol);  //N
// 			println(df.format( 	memb_selfweis[i])+" N");
		}
		constrained = new boolean[nodes.length][];
		for(int i=12;i< nodes.length;i++)
		  constrained[i] = new boolean[] { true, true, true, true, true, true }; // RX, RY, RZ, mx, my, mz;
//		constrained[12] = new boolean[] { true, true, true, true, true, true };
//		constrained[13] = new boolean[] { true, true, true, true, true, true };
//		constrained[15] = new boolean[] { true, true, true, true, true, true };
		
		double[] As = new double[members.length];  //m2, section area
		double[] Es = new double[members.length]; // MPa, elastic modulus, parallel grain
		double[] Gs = new double[members.length];//  MPa, shear modulus, or modulus of rigidity, Logan p259
		double[] Js = new double[members.length]; // m4,  torsion constant, Logan p259
		double[] Iy = new double[members.length]; // m4,  principal moment of inertia
		double[] Iz = new double[members.length]; // m4, principal moment of inertia
		double p4=sec_w * sec_w * sec_w * sec_w;
		for (int i = 0; i < members.length; i++) {
			As[i] = sec_w * sec_w;
			Es[i] = 13E9;  // LI JIan, mu cail ke xue p254, standard for design for timber structures p23
			Gs[i] = 70E6; // G_TR, cross grain, LI JIan, mu cail ke xue p254
			Js[i] = 0.14 * p4; // wiki
			Iy[i] = p4 / 12; // Logan p170
			Iz[i] = p4 / 12;
		}
		fr = new Frame3D(nodes, members, loads, constrained);
		fr.analyzeFrame(As, Es, Gs, Js, Iy, Iz, memb_selfweis);
		postProcess();
	}
	
	private void postProcess(){
		res=new double[members.length][];
		val_max=-1;
		val_min=1e12;
		for (int i = 0; i < members.length; i++) {
			double[] fs = fr.member_forces[i];
			res[i] = new double[] { Math.abs(fs[0]), Math.abs(fs[6]) };
			println("axial " + df.format(fs[0]) + ", " + df.format(fs[6]));
			
//			res[i]=new double[]{Math.sqrt(fs[1] * fs[1] + fs[2] * fs[2]),  Math.sqrt(fs[7] * fs[7] + fs[8] * fs[8])};
//			println("shear "+df.format(Math.sqrt(fs[1] * fs[1] + fs[2] * fs[2])) + ", " + df.format(Math.sqrt(fs[7] * fs[7] + fs[8] * fs[8])));

//			res[i] = new double[] { Math.abs(fs[3]), Math.abs(fs[9]) };
//			println("torsion "+df.format(fs[3]) + ", " + df.format(fs[9]));

//			double m1 = Math.sqrt(fs[4] * fs[4] + fs[5] * fs[5]);
//			double m2 = Math.sqrt(fs[10] * fs[10] + fs[11] * fs[11]);
//			res[i]=new double[]{m1, m2};
//			 println("bend "+df.format(m1) + ", " + df.format(m2));
//			double I = sec_w * sec_w * sec_w * sec_w / 12;
//			println("bend " + df.format(1E-6 * m1 * 0.5 * sec_w / I) + ", " + df.format(1E-6 * m2 * 0.5 * sec_w / I) + " MPa");
			
			for (int j = 0; j < 2; j++) {
				double v = res[i][j];
				if (v > val_max)
					val_max = v;
				if (v < val_min)
					val_min = v;
			}
		}
	}

	public void draw() {
		colorMode(HSB);
		smooth();
		strokeWeight(4);
		background(255);
		noFill();
		stroke(0);
		for (int i = 0; i < nodes.length; i++) {
			if(null!=constrained[i]){
				double[] v= nodes[i];
				pushMatrix();
				translate( (float) (sc*v[0]), (float) (sc*v[1]), (float) (sc*v[2]) );
				box(4);
				popMatrix();
			}
		}
		for (int i = 0; i < members.length; i++) {
			int[] m = members[i];
			double[] pa = nodes[m[0]];
			double[] pb = nodes[m[1]];
			double va = 1- (res[i][0]- val_min) /(val_max- val_min) ;
			double vb = 1- (res[i][1] - val_min) / (val_max- val_min);
			line(pa, (float) va, pb, (float) vb);
		}
	}

	private void line(double[] pa, double[] pb) {
		line((float) (sc*pa[0]), (float) (sc*pa[1]), (float) (sc*pa[2]), (float) (sc*pb[0]), (float) (sc*pb[1]), (float)(sc* pb[2]));
	}
	private void line(double[] pa, float va, double[] pb, float vb) { // va, vb 0-1
		double[] md = M.between(0.5, pa, pb);
		stroke(va * 200, 255, 255);
		line(pa, md);
		stroke(vb * 200, 255, 255);
		line(pb, md);
	}

}
