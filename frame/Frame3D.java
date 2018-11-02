package frame;


public class Frame3D {
	public static final int DIM = 6; // u, v, w, phix, phiy, phiz, each node has 6 D.O.F.
	private double[][] nodes; // x,y,z-coordinates
	private int[][] members; //members[2] ={6,8} means the 3rd member links 7th node to 9th node.
	private double[][][] stiffs, trans; // local stiffness matrix for each member
	private double[][] loads; // nodal: fx, fy, fz, mx, my, mz
	private double[][] eqdis_loads; // equivalent nodal loads given the (uniform) distributed loads
	
	private Double[] memb_selfweis; //self-weight of beam, unit: N or lb.
	private double[] As; // cross-sectional area
	private double[] Es; // Young's modulus
	private double[] Gs; // shear modulus, or modulus of rigidity
	private double[] Js; // torsion constant
	private double[] Iy, Iz; // principal moment of inertia
	private boolean[][] constrained; // constrained[2]={false, true, true} means the 3rd node is 1. constrained along y-axis, and 2. no rotation.
	private Integer[] ni2ai;

	public int lenA; // the number of unconstrained nodal displacements.
	public String[] symbols; // the notations of unconstrained nodal displacements
	// symbols[4]=“2px” means that UA’s 5th element is the 3rd node's rotation around the x-axis.
	private double[][] member_cossin;
	private double[] member_length;
	private double[] member_length2;
	public double[] UA; // (unconstrained) nodal displacement / rotation
	public double[][] member_forces; // the local forces (fx1, fy1, fz1, mx1, my1, mz1, fx2, fy2, fz2, mx2, my2, mz2) of each member.

	public Frame3D(double[][] nodes, int[][] members, double[][] loads, boolean[][] constrained) {
		this.nodes = nodes;
		this.members = members;
		this.loads = loads;
		this.constrained = constrained;
	}

	public void analyzeFrame(double[] As, double[] Es, double[] Gs, double[] Js, double[] Iy, double[] Iz, Double[] memb_selfweis) {
		this.As = As;
		this.Es = Es;
		this.Gs = Gs;
		this.Js = Js;
		this.Iy = Iy;
		this.Iz = Iz;
		this.memb_selfweis = memb_selfweis;
		
		preProcess();
		construct_matrices();
		member_stress();
	}

	private void preProcess() {
		ni2ai = new Integer[DIM* nodes.length];//***** 3: u, v, phi ,   6: u, v, w, phix, phiy, phiz 
		lenA= ni2ac(constrained, ni2ai);
		member_cossin = new double[members.length][];
		member_length = new double[members.length];
		member_length2 = new double[members.length];
		for (int i = 0; i < members.length; i++) {
			double[] pa = nodes[members[i][0]];
			double[] pb = nodes[members[i][1]];
			double x = pb[0] - pa[0];
			double y = pb[1] - pa[1];
			double z = pb[2] - pa[2];
			double len = Math.sqrt(x * x + y * y + z * z);
			member_cossin[i] = new double[] { x / len, y / len, z / len };
			member_length[i] = len;
			member_length2[i] = x * x + y * y + z * z;
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
				symbols[ni2ai[i]]= string6(i);
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
			double[] t = member_cossin[i];
			double[] v = new double[12]; // new double[] { 0, -wf * 0.5, -moment, 0, -wf * 0.5, moment };// 6, right hand rule
			v[2] = v[8] = -wf * 0.5;
			double mx = -wf * t[1] * L / 12.0;  //right hand rule
			double my = wf * t[0] * L / 12.0;  //right hand rule
			v[3] = mx;
			v[4] = my;
			v[9] = -mx;
			v[10] = -my;
			eqdis_loads[i] = v;
			for (int j = 0; j < 6; j++) {
				loads[ida][j] += v[j];
				loads[idb][j] += v[6 + j];
			}
		}
	}

	private static int ni2ac(boolean[][] constrained, Integer[] lista) {
		int active_count = 0;
		for (int i = 0; i < constrained.length; i++) {
			if (null == constrained[i]) {  //proceed DIM steps
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

	private void construct_matrices() {
		double[][] KAA = new double[lenA][lenA];
		for (int i = 0; i <members.length; i++) { // members.length
			int ida = members[i][0];
			int idb = members[i][1];
			if (ida > idb)
				throw new RuntimeException();
			double[][] v= trans[i];
			double[][] mat = M.mul(M.transpose(v), M.mul(stiffs[i], v)); // 12 * 12
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
		UA = M.solve_square_Axb(KAA, FA);  //use Jama to solve. One might replace this method using other linear solvers.
	}

	private double[][] stiff(int i) { // member id, local stiffness, Hutton 278
		double L = member_length[i];
		double L2 = member_length2[i];
		double L3 = L * L2;
		double A = As[i];
		double E = Es[i];
		double iz = Iz[i];
		double iy = Iy[i];
		double G = Gs[i];
		double J = Js[i];
		double[][] ker = new double[2 * DIM][2 * DIM]; // 12*12
		
		ker[0][0] = ker[DIM][DIM] = A * E / L; // ******************************axial
		ker[0][DIM] = ker[DIM][0] = -A * E / L;

		ker[1][1] = ker[1 + DIM][1 + DIM] = 12 * E * iz / L3; // ******************** (frame) shear bend
		ker[1][1 + DIM] = ker[1 + DIM][1] = -12 * E * iz / L3;
		ker[1][5] = ker[1][5 + DIM] = ker[5][1] = ker[5 + DIM][1] = 6 * E * iz / L2;
		ker[1 + DIM][5] = ker[1 + DIM][5 + DIM] = ker[5][1 + DIM] = ker[5 + DIM][1 + DIM] = -6 * E * iz / L2;
		ker[5][5] = ker[5 + DIM][5 + DIM] = 4 * E * iz / L;
		ker[5][5 + DIM] = ker[5 + DIM][5] = 2 * E * iz / L;

		ker[2][2] = ker[2 + DIM][2 + DIM] = 12 * E * iy / L3;// ******************** (grid) torsion shear bend
		ker[2][2 + DIM] = ker[2 + DIM][2] = -12 * E * iy / L3;
		ker[2][4] = ker[2][4 + DIM] = ker[4][2] = ker[4 + DIM][2] = -6 * E * iy / L2;//
		ker[2 + DIM][4] = ker[2 + DIM][4 + DIM] = ker[4][2 + DIM] = ker[4 + DIM][2 + DIM] = 6 * E * iy / L2;
		ker[4][4] = ker[4 + DIM][4 + DIM] = 4 * E * iy / L;
		ker[4][4 + DIM] = ker[4 + DIM][4] = 2 * E * iy / L;
		ker[3][3] = ker[3 + DIM][3 + DIM] = G * J / L;
		ker[3 + DIM][3] = ker[3][3 + DIM] = -G * J / L;
		return ker;
	}

	private double[][] tranform(int i) { // member id
		double[] xp = member_cossin[i]; // x/l, y/l, z/l
		double[] yp;
		if (M.concide(new double[] { 0, 0, 1 }, xp)) {
			//println("concide");
			yp = new double[] { 0, 1, 0 };
		} else {
			yp = M.cross(new double[] { 0, 0, 1 }, xp);
			M._normalize(yp);
		}
		double[] zp = M.cross(xp, yp);
		double[][] T = { xp, yp, zp };
		double[][] v = new double[2 * DIM][2 * DIM]; // 12*12
		for (int k = 0; k < 4; k++) {
			for (int m = 0; m < 3; m++)
				for (int n = 0; n < 3; n++)
					v[m + 3 * k][n + 3 * k] = T[m][n];
		}
		return v;
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
			double[] u = new double[2 * DIM];//12
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
				M._sub(member_forces[i], M.mul(trans[i], eqdis_loads[i]));//distributed load
		}
	}

	private static String string6(int i) {
		String t;
		if (0 == i % DIM)
			t = "u";
		else if (1 == i % DIM)
			t = "v";
		else if (2 == i % DIM)
			t = "w";
		else if (3 == i % DIM)
			t = "px";
		else if (4 == i % DIM)
			t = "py";
		else 
			t = "pz";
		t =(i / DIM)+t;
		return t;
	}
}
