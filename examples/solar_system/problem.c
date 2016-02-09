/**
 * Solar System
 *
 * This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 Log:
 Jan 14: 
>>we are assuming velocity in calculation is Barycentric velocity and the force calculation(the e.o.m) is 
using the veloctiy relativte to the Sun, i.e. v_i - v_s.
>>and in the H_pn, we are certain that those v_tilde  are velocity relative to the Sun.
>> Question:
 >>One thing we can be sure is that, for calculating the force, we have to assume the velocity in computation 
 >>is Barycentric and use v_i-v_S. We are trying to get the whole energy of the system, so think of H_pn as just 
 >>the relativistic correction, so we break the total energy into two parts, first is the classical energy where
 >>it includes all planets and Sun's kinetic energy and plus the potential energy between planets and the Sun.
 >>Second one is H_pn, where the v_tilde in the expression should be just velocity relative to the Sun, since we 
 >>we derived the same e.o.m from this H_pn as given by Benitez&Gallardo which has v_tilde as relative velocity.

 Jan 21:
 >> Doing only p dependent H: p^4
 Feb 3:
 >> remove or add the back reaction term in force calculation seems not affecting the error too much, and also
 + or - sign doesn't seem like affecting the result.
 >> even if I multiplied 1e3 on pi.ax += ...1.e3, or just simply take pi.ax += 1, it didn't affect the result.
 >> the relativistic correction printed to be around 1e-11 ~ 1e-13.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double ss_pos[10][3] = 
{
	{3.256101656448802E-03  , -1.951205394420489E-04 , -1.478264728548705E-04},  
	// {-1.927589645545195E-01 , 2.588788361485397E-01  , 3.900432597062033E-02 }, 
	// {-5.976537074581466E-01 , 3.918678996109574E-01  , 3.990356741282203E-02 }, 
	// // {-7.986189029000561E-01 , -6.086873314992410E-01 , -1.250824315650566E-04}, 
	// // {7.897942807177173E-01  , 1.266671734964037E+00  , 7.092292179885432E-03 }, 
	{-4.314503046344270E+00, 3.168094294126697E+00, 8.331048545353310E-02 }, 
	// {-4.882304833383455E+00 , -8.689263067189865E+00 , 3.453930436208210E-01 }, 
	// {1.917757033372740E+01  , 5.671738750949031E+00  , -2.273858614425555E-01},  
	// {2.767031517959636E+01  , -1.150331645280942E+01 , -4.008018419157927E-01},  
	// {7.765250227278298E+00  , -3.190996242617413E+01 , 1.168394015703735E+00 }, 

};
double ss_vel[10][3] = 
{
	{3.039963463108432E-06 ,  6.030576499910942E-06 ,  -7.992931269075703E-08}, 
	// {-2.811550184725887E-02,  -1.586532995282261E-02,  1.282829413699522E-03 }, 
	// {-1.113090630745269E-02,  -1.703310700277280E-02,  4.089082927733997E-04 },
	// // {1.012305635253317E-02 ,  -1.376389620972473E-02,  3.482505080431706E-07 }, 
	// // {-1.135279609707971E-02,  8.579013475676980E-03 ,  4.582774369441005E-04 }, 
	{-4.555986691913995E-03,  -5.727124269621595E-03,  1.257262404884127E-04 }, 
	// {4.559352462922572E-03 ,  -2.748632232963112E-03,  -1.337915989241807E-04}, 
	// {-1.144087185031310E-03,  3.588282323722787E-03 ,  2.829006644043203E-05 }, 
	// {1.183702780101068E-03 ,  2.917115980784960E-03 ,  -8.714411604869349E-05}, 
	// {3.112825364672655E-03 ,  1.004673400082409E-04 ,  -9.111652976208292E-04},
};

double ss_mass[10] =
{
	1.988544e30,
	// 3.302e23,
	// 48.685e23,
	// // 6.0477246e24,
	// // 6.4185e23,
	1898.13e24,
	// 5.68319e26,
	// 86.8103e24,
	// 102.41e24,
	// 1.4639248e+22,
};

void heartbeat(struct reb_simulation* r);
double energy(struct reb_simulation* r);
void additional_forces(struct reb_simulation* r);

double e_init;
double tmax;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->dt 			= 4.;				// in days
	tmax			= 365.25 * 1.e4/2. ;//7.3e10;			// 200 Myr
	r->G			= 1.4880826e-34;		// in AU^3 / kg / day^2.
	//r->ri_whfast.safe_mode 	= 0;		// Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
	//r->ri_whfast.corrector 	= 11;		// 11th order symplectic corrector
	// r->integrator		= REB_INTEGRATOR_WHFAST;
	r->heartbeat		= heartbeat;
	r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
	r->integrator		= REB_INTEGRATOR_IAS15;		// Alternative non-symplectic integrator

	r->additional_forces 	= additional_forces;
	r->force_is_velocity_dependent = 1;

	// Initial conditions
	for (int i=0;i<2;i++){
		struct reb_particle p = {0};
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.m  = ss_mass[i];
		reb_add(r, p); 
	}
	reb_move_to_com(r);
	e_init = energy(r);
	system("rm -f energy.txt");
	reb_integrate(r, tmax);
}



// Simple(non-complete) version // 10e-14
void additional_forces(struct reb_simulation* r){
	struct reb_particle* particles = r->particles;
	const int N = r->N;
	double G = r->G;
	double C = 3.e8*24.*60.*60./(149.6e6 * 1000);// * 1.e3; // AU / day
	// struct reb_particle com = reb_get_com(r)
	for (int i=1; i<N; i++){ 
		// compute only relative to the Sun, in this case it's pulsar A
		struct reb_particle sun = particles[0];
		struct reb_particle pi = particles[i];
		double dx = particles[i].x-sun.x; // the distance is relative to sun
		double dy = particles[i].y-sun.y;
		double dz = particles[i].z-sun.z;
		double r2 = dx*dx + dy*dy + dz*dz;
		double r = sqrt(r2);
		double dvx = particles[i].vx;//-sun.vx; // note: the velocity is now chnaged to veloctiy relative to sun
		double dvy = particles[i].vy;//-sun.vy;
		double dvz = particles[i].vz;//-sun.vz;
		// double diff = pi.ax - (-G*sun.m/(r2*r)*dx);
		// printf(">>>>>%.15e >>>>>%.15e >>>>>%.15e\n", diff, pi.ax, G*sun.m/(r2*r)*dx);
		
		// Benitez and Gallardo 2008
		double alpha = G*sun.m/(r*r*r*C*C);
		double v_mag = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
		double beta = v_mag*v_mag/2.;
		double gamma = (dvx*dx + dvy*dy + dvz*dz);

		particles[i].ax += alpha*(beta*dx + gamma*dvx);
		particles[i].ay += alpha*(beta*dy + gamma*dvy);
		particles[i].az += alpha*(beta*dz + gamma*dvz);
		particles[0].ax -= particles[i].m/particles[0].m*alpha*(beta*dx + gamma*dvx);
		particles[0].ay -= particles[i].m/particles[0].m*alpha*(beta*dy + gamma*dvy);
		particles[0].az -= particles[i].m/particles[0].m*alpha*(beta*dz + gamma*dvz);
		// double ppx = 0.;
		// double ppy = 0.;
		// double ppz = 0.;
		// double px = pi.m*vx;
		// double py = pi.m*vy;
		// double pz = pi.m*vz;

		// double pp2 = 0.;
		// double p2 = px*px + py*py + pz*pz;
		// while (fabs(p2 - pp2)>1.e-15){
		// 	pp2 =p2;
		// 	ppx = px; ppy = py; ppz = pz;
		// 	px = pi.m*vx + pp2*ppx/(2.*pi.m*pi.m*C*C);
		// 	py = pi.m*vy + pp2*ppy/(2.*pi.m*pi.m*C*C);
		// 	pz = pi.m*vz + pp2*ppz/(2.*pi.m*pi.m*C*C);
		// 	p2 = px*px + py*py + pz*pz;
		// }
		// double pdot = px*pi.m*pi.ax + py*pi.m*pi.ay + pz*pi.m*pi.az;
		// pi.ax += (2.*pdot*px + p2*pi.m*pi.ax)/(2.*pi.m*pi.m*pi.m*C*C);
		// pi.ay += (2.*pdot*py + p2*pi.m*pi.ay)/(2.*pi.m*pi.m*pi.m*C*C);
		// pi.az += (2.*pdot*pz + p2*pi.m*pi.az)/(2.*pi.m*pi.m*pi.m*C*C);
		// // // printf(">>>>>>>>>>>>>>>>>>>>>>%.15e\n", (2.*pdot*px + p2*pi.m*pi.ax)/(2.*pi.m*pi.m*pi.m*C*C));
		// sun.ax -= pi.m/sun.m*(2.*pdot*px + p2*pi.m*pi.ax)/(2.*pi.m*pi.m*pi.m*C*C);
		// sun.ay -= pi.m/sun.m*(2.*pdot*py + p2*pi.m*pi.ay)/(2.*pi.m*pi.m*pi.m*C*C);
		// sun.az -= pi.m/sun.m*(2.*pdot*pz + p2*pi.m*pi.az)/(2.*pi.m*pi.m*pi.m*C*C);

	}
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Simplified PN potential energy //10e-14
double energy(struct reb_simulation* r){
	struct reb_particle* particles = r->particles;
	int N = r->N;
	double G = r->G;
	double C = 3.e8*24.*60.*60./(149.6e6 * 1000);// * 1.e3; // AU / day
	double K2 = G*ss_mass[0];

	double e_kin = 0.;
	double e_pot = 0.;
	double e_pn  = 0.;
	struct reb_particle sun = particles[0];
	for (int i=0;i<N;i++){
		struct reb_particle pi = particles[i];
		if (i != 0){
			double vx = pi.vx;// - sun.vx;
			double vy = pi.vy;// - sun.vy;
			double vz = pi.vz;// - sun.vz;
			double v2 = vx*vx + vy*vy + vz*vz;
			
			double A = 1. - v2/(2.*C*C);
			// double p = sqrt(v2)/A * pi.m;
			double B = sqrt(v2)/A;
			double v_tilde2 = B*B;
			// double ppx = 0.;
			// double ppy = 0.;
			// double ppz = 0.;
			// double px = pi.m*vx;
			// double py = pi.m*vy;
			// double pz = pi.m*vz;

			// double pp2 = 0.;
			// double p2 = px*px + py*py + pz*pz;
			// while (fabs(p2 - pp2)>1.e-15){ // comment out the while loop if classical velocity needed
			// 	pp2 =p2;
			// 	ppx = px; ppy = py; ppz = pz;
			// 	px = pi.m*vx + pp2*ppx/(2.*pi.m*pi.m*C*C);
			// 	py = pi.m*vy + pp2*ppy/(2.*pi.m*pi.m*C*C);
			// 	pz = pi.m*vz + pp2*ppz/(2.*pi.m*pi.m*C*C);
			// 	p2 = px*px + py*py + pz*pz;
			// }
			// e_kin += 0.5*p2/pi.m;
			// e_pn += (1./(C*C)) * (-p2*p2/(8.*pi.m*pi.m*pi.m));
			// e_kin += 0.5*pi.m*v_tilde2;
			e_kin += 0.5*v_tilde2;
			// e_pn -= v_tilde2*v_tilde2*pi.m/(8.*C*C);
			e_pn -= v_tilde2*v_tilde2/(8.*C*C);
		}		
		else{
			double sun_v2 = sun.vx*sun.vx + sun.vy*sun.vy + sun.vz*sun.vz;
			// e_kin += 0.5 * sun.m * sun_v2;
			e_kin += 0.5 * (sun.m/ss_mass[1]) * sun_v2;
		}
		for (int j=i+1; j<N; j++){
			struct reb_particle pj = particles[j];
			double dx = pi.x - pj.x;
			double dy = pi.y - pj.y;
			double dz = pi.z - pj.z;	
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			
			// e_pot -= G*pi.m*pj.m/r;
			e_pot -= G*pi.m/r;
		}
	}
	return e_kin + e_pot + e_pn;
}


// classic Newtonian energy
// double energy(struct reb_simulation* const r){
// 	struct reb_particle* const particles = r->particles;
// 	const int N = r->N;
// 	double G = r->G;
// 	double e_kin = 0.;
// 	double e_pot = 0.;
// 	for (int i=0;i<N;i++){
// 		struct reb_particle pi = particles[i];
// 		e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
// 		for (int j=i+1;j<N;j++){
// 			struct reb_particle pj = particles[j];
// 			double dx = pi.x - pj.x;
// 			double dy = pi.y - pj.y;
// 			double dz = pi.z - pj.z;
// 			e_pot -= G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
// 		}
// 	}
// 	return e_kin +e_pot;
// }

double next_output = 1.;
void heartbeat(struct reb_simulation* r){
	if (r->t > next_output){
		next_output *= 1.02;
		struct reb_particle com = reb_get_com(r);
		double x = com.x;
		double y = com.y;
		double z = com.z;
		double r2 = x*x + y*y + z*z;
		// double r = sqrt(r2);
		// struct reb_particle com = reb_get_com(r);
		// printf("%e, %e, %e\n", com.x, com.y, com.z);
		//reb_output_timing(r, tmax);
		reb_integrator_synchronize(r);
		FILE* f = fopen("energy.txt","a");
		double e = energy(r);
		fprintf(f, "%e %.16e\n", r->t, fabs((e-e_init)/e_init));
		// fprintf(f, "%e %.16e\n", r->t, sqrt(r2));
		fclose(f);
	}
}

