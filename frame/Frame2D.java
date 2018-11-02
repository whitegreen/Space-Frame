package frame;

//Hao Hua, Southeast University, whitegreen@163.com

public class Frame2D {
	public static final  int DIM=3;  // u, v, pi, each node has 3 D.O.F. 
	private double[][] nodes;   //x,y-coordinates
	private int[][] members ; // members[2] ={6,8} means the 3rd member links 7th node to 9th node.
	private double[][][] stiffs, trans;  //local stiffness matrix for each member
	private double[][] loads;  //nodal, fx, fy, bending
	private double[][] eqdis_loads;  // equivalent nodal loads given the  (uniform) distributed loads
	private Double[] memb_selfweis;   //self-weight of beam, unit: N or lb.
	
	private double[] As; // cross-sectional area 
	private double[] Es; // Young's modulus
	private double[] Is;  //principal moment of inertia
	private boolean[][] constrained; // constrained[2]={false, true, true} means the 3rd node is 1. constrained along y-axis, and 2. no rotation.
	private Integer[] ni2ai;// 1D node index (DIM*) to active index

	public int lenA;   //the number of unconstrained nodal displacements.
	public String[] symbols;  //the notations of unconstrained nodal displacements
	//symbols[4]=“2p” means that UA’s 5th element is the 3rd node's rotation.
	private double[][] member_cossin;
	private double[] member_length;
	public double[] UA; // (unconstrained) nodal displacement / rotation
	public double[][] member_forces;  //the local forces (fx1, fy1, m1, fx2,fy2, m2) of each member.
	
	public Frame2D(double[][] nodes, int[][] members, double[][] loads, boolean[][] constrained, double[] As, double[] Es, double[] Is, Double[] memb_selfweis) {
		this.nodes = nodes;
		this.members = members;
		this.loads = loads;
		this.constrained = constrained;
		this.As = As;
		this.Es = Es;
		this.Is = Is;
		this.memb_selfweis = memb_selfweis;

		preProcess();
		construct_matrices();
		member_stress();
	}

	private void preProcess() {
		ni2ai = new Integer[DIM * nodes.length];// ***** 3: u, v, phi
		lenA = ni2ac(constrained, ni2ai);
		member_cossin = new double[members.length][];
		member_length = new double[members.length];
		for (int i = 0; i < members.length; i++) {
			double[] pa = nodes[members[i][0]];
			double[] pb = nodes[members[i][1]];
			double x = pb[0] - pa[0];
			double y = pb[1] - pa[1];
			double len = Math.sqrt(x * x + y * y);
			member_cossin[i] = new double[] { x / len, y / len };
			member_length[i] = len;
		}
		stiffs = new double[members.length][][];
		trans = new double[members.length][][];
		for (int i = 0; i < members.length; i++) {
			stiffs[i] = stiff(i);
			trans[i] = tranform(i);
		}
		if (null != memb_selfweis)
			setDistributedLoads();
		
		symbols=new String[lenA];
		for (int i = 0; i < ni2ai.length; i++) {
			if (null != ni2ai[i]) 
				symbols[ni2ai[i]]= string3(i);
		}
	}

	private void setDistributedLoads() {
		eqdis_loads = new double[members.length][];
		for (int i = 0; i < members.length; i++) {
			Double wf = memb_selfweis[i];
			if (null == wf)
				continue;
			double L = member_length[i];
			int ida = members[i][0];
			int idb = members[i][1];
			if (ida > idb)
				throw new RuntimeException();
			double moment = wf * L * member_cossin[i][0] / 12.0;
			double[] v = new double[] { 0, -wf * 0.5, -moment, 0, -wf * 0.5, moment };// 6
			eqdis_loads[i] = v;
			for (int j = 0; j < 3; j++) {
				loads[ida][j] += v[j];
				loads[idb][j] += v[3 + j];
			}
		}
	}

	private static int ni2ac(boolean[][] constrained, Integer[] lista) {
		int active_count = 0;
		for (int i = 0; i < constrained.length; i++) {
			if (null == constrained[i]) { // proceed DIM steps
				for (int k = 0; k < DIM; k++)
					lista[DIM * i + k] = active_count + k;
				active_count += DIM;
				continue;
			}
			for (int k = 0; k < DIM; k++) {
				if (!constrained[i][k]) {
					lista[DIM * i + k] = active_count;
					active_count++;
				}
			}
		}
		return active_count;
	}

	private double[][] stiff(int i) { // member id, local stiffness, Hutton p116, Logan p238
		double L = member_length[i];
		double L2 = L * L;
		double ae = As[i] * Es[i] / L;
		double ei = Es[i] * Is[i] / L;
		double[][] m = new double[6][6];
		m[0][0] = m[3][3] = ae;
		m[0][3] = m[3][0] = -ae;
		m[1][1] = m[4][4] = 12 * ei / L2;
		m[1][4] = m[4][1] = -12 * ei / L2;
		m[1][2] = m[1][5] = m[2][1] = m[5][1] = 6 * ei / L;
		m[4][2] = m[2][4] = m[4][5] = m[5][4] = -6 * ei / L;
		m[2][2] = m[5][5] = 4 * ei;
		m[2][5] = m[5][2] = 2 * ei;
		return m;
	}

	private double[][] tranform(int i) { // member id
		double[] t = member_cossin[i]; // cos, sin
		double[][] v = new double[6][6];
		v[0][0] = v[1][1] = v[3][3] = v[4][4] = t[0]; // c
		v[2][2] = v[5][5] = 1;
		v[0][1] = v[3][4] = t[1]; // s
		v[1][0] = v[4][3] = -t[1];
		return v;
	}

	private void construct_matrices() {
		double[][] KAA = new double[lenA][lenA];
		for (int i = 0; i < members.length; i++) { // members.length
			int ida = members[i][0];
			int idb = members[i][1];
			if (ida > idb)
				throw new RuntimeException();
			double[][] v = trans[i];
			double[][] mat = M.mul(M.transpose(v), M.mul(stiffs[i], v));
			for (int m = 0; m < DIM; m++) {
				for (int n = 0; n < DIM; n++) {
					assemble(DIM * ida + m, DIM * ida + n, ni2ai, mat[m][n], KAA);
					assemble(DIM * ida + m, DIM * idb + n, ni2ai, mat[m][DIM + n], KAA);
					assemble(DIM * idb + m, DIM * ida + n, ni2ai, mat[DIM + m][n], KAA);
					assemble(DIM * idb + m, DIM * idb + n, ni2ai, mat[DIM + m][DIM + n], KAA);
				}
			}
		}
		double[] FA = new double[lenA];
		for (int i = 0; i < loads.length; i++) {
			if (null == loads[i])
				continue;
			for (int k = 0; k < DIM; k++) {
				Integer id = ni2ai[i * DIM + k];
				if (null != id)
					FA[id] = loads[i][k];
			}
		}
		UA = M.solve_square_Axb(KAA, FA); //use Jama to solve. One might replace this method using other linear solvers.
	}

	private static void assemble(int a, int b, Integer[] ni2ai, double v, double[][] KAA) {
		Integer ia = ni2ai[a];
		Integer ib = ni2ai[b];
		if (null != ia && null != ib) 
			KAA[ia][ib] += v;
	}

	private void member_stress() {
		member_forces = new double[members.length][];
		for (int i = 0; i < members.length; i++) {
			int ida = members[i][0];
			int idb = members[i][1];
			if (ida > idb)
				throw new RuntimeException();
			double[] u = new double[2 * DIM];// 6
			for (int j = 0; j < DIM; j++) {
				Integer ia = ni2ai[DIM * ida + j];
				if (null != ia)
					u[j] = UA[ia];
				Integer ib = ni2ai[DIM * idb + j];
				if (null != ib)
					u[DIM + j] = UA[ib];
			}
			member_forces[i] = M.mul(stiffs[i], M.mul(trans[i], u));
			if (null != memb_selfweis && null != memb_selfweis[i])
				M._sub(member_forces[i], M.mul(trans[i], eqdis_loads[i]));
		}
	}

	private static String string3(int i) {
		String t;
		if (0 == i % DIM)
			t = "u";
		else if (1 == i % DIM)
			t = "v";
		else
			t = "p";
		t =(i / DIM)+t;
		return t;
	}
	
	
}
