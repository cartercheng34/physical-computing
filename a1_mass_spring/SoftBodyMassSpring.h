//#####################################################################
// Mass-spring deformable model
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################

#ifndef __SoftBodyMassSpring_h__
#define __SoftBodyMassSpring_h__
#include "Common.h"
#include "Particles.h"
#include <typeinfo>

template<int d> class SoftBodyMassSpring
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using MatrixD=Matrix<real,d>;
public:
	////Spring parameters
	Particles<d> particles;
	Array<Vector2i> springs;
	Array<real> rest_length;
	Array<real> ks;
	Array<real> kd;

	////Boundary nodes
	Hashtable<int,VectorD> boundary_nodes;

	////Body force
	VectorD g=VectorD::Unit(1)*(real)-1.;
	
	enum class TimeIntegration{ExplicitEuler,ImplicitEuler} time_integration=TimeIntegration::ImplicitEuler;

	////Implicit time integration
	SparseMatrixT K, Mass;
	VectorX u,b;

	virtual void Initialize()
	{
		////Initialize default spring parameters for standard tests
		real ks_0=(real)1,kd_0=(real)1;
		switch(time_integration){
		case TimeIntegration::ExplicitEuler:{
			ks_0=(real)5e2;
			kd_0=(real)1e1;
		}break;
		case TimeIntegration::ImplicitEuler:{
			ks_0=(real)2e5;
			kd_0=(real)1e1;			
		}break;}

		////Allocate arrays for springs and parameters
		rest_length.resize(springs.size());
		for(int i=0;i<(int)springs.size();i++){const Vector2i& s=springs[i];
			rest_length[i]=(particles.X(s[0])-particles.X(s[1])).norm();}
		ks.resize(springs.size(),ks_0);
		kd.resize(springs.size(),kd_0);

		////Allocate sparse matrix if using implicit time integration 
		////This function needs to be called for only once since the mesh doesn't change during the simulation)
		if(time_integration==TimeIntegration::ImplicitEuler)
			Initialize_Implicit_K_And_b();
	}

	virtual void Advance(const real dt)
	{
		switch(time_integration){
		case TimeIntegration::ExplicitEuler:
			Advance_Explicit_Euler(dt);break;
		case TimeIntegration::ImplicitEuler:
			Advance_Implicit_Euler(dt);break;}
	}
	
	////Set boundary nodes
	void Set_Boundary_Node(const int p,const VectorD v=VectorD::Zero()){boundary_nodes[p]=v;}
	
	bool Is_Boundary_Node(const int p){return boundary_nodes.find(p)!=boundary_nodes.end();}
	
	void Enforce_Boundary_Conditions()
	{
		for(auto p:boundary_nodes){
			int idx=p.first;					////get boundary particle index
			const VectorD& v=p.second;			////get boundary particle velocity
			particles.V(idx)=v;					////set boundary particle velocity
			particles.F(idx)=VectorD::Zero();}	////clear boundary particle force
	}

	//////////////////////////////////////////////////////////////////////////
	////P1 TASK: explicit Euler integration and spring force calculation

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): explicit Euler time integration 
	void Advance_Explicit_Euler(const real dt)
	{
		Particle_Force_Accumulation();

		////Update particle velocity and position
		/* Your implementation start */
		for(int i = 0 ;i < particles.Size() ; i++){
			particles.V(i) += particles.F(i) / particles.M(i) * dt;
			particles.X(i) += particles.V(i)*dt;
		}
		/* Your implementation end */
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute spring force f_ij=f_s+f_d 
	VectorD Spring_Force_Calculation(const int idx)
	{
		/* Your implementation start */
		
		int idx_i = springs[idx][0];
		int idx_j = springs[idx][1];
		auto unit_v = (particles.X(idx_j) -  particles.X(idx_i)).normalized();
		auto fs_i = ks[idx] * ( (particles.X(idx_j) -  particles.X(idx_i)).norm() -  rest_length[idx]) * unit_v ;
		auto fd_i = kd[idx] * (particles.V(idx_j) -  particles.V(idx_i)).dot(unit_v) * unit_v;
		auto f_i = fs_i + fd_i;
		return f_i; 
		
			
		
		/* Your implementation end */

		//return VectorD::Zero();	////replace this line with your implementation
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): accumulate spring forces to particles
	void Particle_Force_Accumulation()
	{
		////Clear forces on particles
		for(int i=0;i<particles.Size();i++){particles.F(i)=VectorD::Zero();}

		////Accumulate body forces
		for(int i=0;i<particles.Size();i++){
			particles.F(i) += particles.M(i)*g;
		}

		////Accumulate spring forces
		/* Your implementation start */
		for (int i = 0 ;i < springs.size() ; i++){
			auto f_i = Spring_Force_Calculation(i);
			auto f_j = -f_i;
			int idx_i = springs[i][0];
			int idx_j = springs[i][1];
			particles.F(idx_i) += f_i;
			particles.F(idx_j) += f_j;
		}
		/* Your implementation end */

		////Enforce boundary conditions
		Enforce_Boundary_Conditions();
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 1: initialize the matrix structure 
	void Initialize_Implicit_K_And_b()
	{
		int n=d*particles.Size();
		Mass.resize(n, n);
		K.resize(n,n);u.resize(n);u.fill((real)0);b.resize(n);b.fill((real)0);
		Array<TripletT> elements;
		for(int s=0;s<(int)springs.size();s++){int i=springs[s][0];int j=springs[s][1];
			Add_Block_Triplet_Helper(i,i,elements);
			Add_Block_Triplet_Helper(i,j,elements);
			Add_Block_Triplet_Helper(j,i,elements);
			Add_Block_Triplet_Helper(j,j,elements);}
		K.setFromTriplets(elements.begin(),elements.end());
		K.makeCompressed();	
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2: fill nonzero elements in K
	void Update_Implicit_K_And_b(const real dt)
	{
		////Clear K and b
		K.setZero();
		b.fill((real)0);
		Mass.setZero();
		

		VectorXd total_f(d*particles.Size());
		VectorXd total_v(d*particles.Size());
		
		// init Mass matrix as sparse matrix
		for (int i = 0 ; i < particles.Size() ; i++){
			// for (int j = 0 ; j < d ; j++){
			// 	total_v[i*d + j] = particles.V(i)[j];
			// 	total_f[i*d + j] = particles.F(i)[j];
			// }
			Set_Block(total_f, i, particles.F(i));
			Set_Block(total_v, i, particles.V(i));
			
			auto tmp_M1 = MatrixD::Identity() * particles.M(i);
			SparseFunc::Add_Block<d,MatrixD>(Mass,i,i,tmp_M1);
		}
		
		
		/* Your implementation start */
		for (int i = 0; i < springs.size() ; i++){
			MatrixD damping_m;
			MatrixD stiff_m;
			MatrixXd Mass(6,6);
			
			auto tmp_M1 = MatrixD::Identity() * particles.M(springs[i][0]);
            auto tmp_M2 = MatrixD::Identity() * particles.M(springs[i][1]);
			
			Mass << tmp_M1, MatrixD::Zero(), MatrixD::Zero(), tmp_M2 ;
            
            
			Compute_Kd_Block(i, damping_m);
			Compute_Ks_Block(i, stiff_m);

			Add_Block_Helper(K, springs[i][0], springs[i][1], -damping_m * dt - stiff_m * dt *dt);
			
			// Add_Mass_Helper(K, springs[i][0], tmp_M1);
			// Add_Mass_Helper(K, springs[i][1], tmp_M2);


			VectorXd tmp_v(6);
			tmp_v << particles.V(springs[i][0]), particles.V(springs[i][1]);
			
			// VectorXd tmp_b = Mass * tmp_v;

			// VectorXd tmp_f(6);
			// tmp_f << particles.F(springs[i][0]), particles.F(springs[i][1]);
			// tmp_b += tmp_f * dt;

			MatrixXd tmp_d(6,6);
			tmp_d << damping_m, -damping_m, -damping_m, damping_m;
			VectorXd tmp_b = -tmp_d * tmp_v * dt;

			Add_Block(b, springs[i][0], tmp_b.head(3));
			Add_Block(b, springs[i][1], tmp_b.tail(3));
			
		}

		K += Mass;
		b += Mass * total_v;
		b += total_f * dt;

		/* Your implementation end */
	}

	//////////////////////////////////////////////////////////////////////////
	////P2 TASK: Implicit Euler time integration

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2.1: compute spring force derivative
	void Compute_Ks_Block(const int s,MatrixD& Ks)
	{
		/* Your implementation start */
		int idx_i = springs[s][0];
		int idx_j = springs[s][1];
		auto dir = particles.X(idx_j) -  particles.X(idx_i);
		auto unit_v = (particles.X(idx_j) -  particles.X(idx_i)).normalized();
		auto distance = (particles.X(idx_j) -  particles.X(idx_i)).norm();
		auto K_ii = ks[s] * ((rest_length[s]/distance - 1) * Matrix3d::Identity() - rest_length[s]/pow(distance, 3) * dir * dir.transpose());
		
		Ks = K_ii;
		/* Your implementation end */
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2.2: compute damping force derivative
	void Compute_Kd_Block(const int s,MatrixD& Kd)
	{
		/* Your implementation start */
		int idx_i = springs[s][0];
		int idx_j = springs[s][1];
		auto dir = particles.X(idx_j) -  particles.X(idx_i);
		auto unit_v = (particles.X(idx_j) -  particles.X(idx_i)).normalized();
		auto distance = (particles.X(idx_j) -  particles.X(idx_i)).norm();
		auto D_ii = -kd[s] * dir * dir.transpose()/ pow(distance, 2);
		
		Kd = D_ii;
		
		
		/* Your implementation end */
	}

	////Implicit Euler time integration
	void Advance_Implicit_Euler(const real dt)
	{
		Particle_Force_Accumulation();
		Update_Implicit_K_And_b(dt);

		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<d;j++) u[i*d+j]=particles.V(i)[j];}	////set initial guess to be the velocity from the last time step

		SparseSolver::CG(K,u,b);	////solve Ku=b using Conjugate Gradient

		for(int i=0;i<particles.Size();i++){
			VectorD v;for(int j=0;j<d;j++) v[j]=u[i*d+j];
			particles.V(i)=v;
			particles.X(i)+=particles.V(i)*dt;}
	}

protected:
	////Add block nonzeros to sparse matrix elements (for initialization)
	void Add_Block_Triplet_Helper(const int i,const int j,Array<TripletT>& elements)
	{for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++)elements.push_back(TripletT(i*d+ii,j*d+jj,(real)0));}

	////Add block Ks to K_ij
	void Add_Block_Helper(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
		SparseFunc::Add_Block<d,MatrixD>(K,j,j,Ks);
		if(!Is_Boundary_Node(i)&&!Is_Boundary_Node(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,i,j,-Ks);
			SparseFunc::Add_Block<d,MatrixD>(K,j,i,-Ks);}
	}

	void Add_Mass_Helper(SparseMatrixT& K,const int i,const MatrixD& Ks)
	{
		if(!Is_Boundary_Node(i))
			SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
	}

	////Set block values on a vector
	void Set_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]=bi[ii];}

	////Add block values to a vector
	void Add_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]+=bi[ii];}
};

#endif
