#include <cmath>
#include <iostream>
#include <limits>
#include <random>

#if __has_include("pcg/pcg_random.hpp")
#include "pcg/pcg_random.hpp"
#define USE_PCG
#elif __has_include("pcg_random.hpp")
#define USE_PCG
#endif

#include "cl_options.h"

struct Position{
	double x,y,z;
	Position():x(0),y(0),z(0){}
	Position(double x, double y, double z):x(x),y(y),z(z){}
	Position(const Position& p):x(p.x),y(p.y),z(p.z){}
	
	double magnitude() const{ return sqrt(x*x+y*y+z*z); }
	double mag2() const{ return (x*x+y*y+z*z); }
	
	Position& operator+=(const Position& p){
		x+=p.x;
		y+=p.y;
		z+=p.z;
		return *this;
	}
	Position& operator-=(const Position& p){
		x-=p.x;
		y-=p.y;
		z-=p.z;
		return *this;
	}
	
	Position operator+(const Position& p) const{
		return Position(*this)+=p;
	}
	Position operator-(const Position& p) const{
		return Position(*this)-=p;
	}
	
	Position& operator*=(double a){
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	Position& operator/=(double a){
		x/=a;
		y/=a;
		z/=a;
		return *this;
	}
	
	Position operator*(double a) const{
		return Position(*this)*=a;
	}
	Position operator/(double a) const{
		return Position(*this)/=a;
	}
	
	Position& normed(){
		return *this/=magnitude();
	}
	Position unit() const{
		return Position(*this).normed();
	}
	
	double operator*(const Position& p) const{
		return x*p.x+y*p.y+z*p.z;
	}
	double dot(const Position& p) const{
		return x*p.x+y*p.y+z*p.z;
	}
};

Position operator*(double a, const Position& p){
	return Position(p)*=a;
}
double dot(const Position& p1, const Position& p2){
	return p1.dot(p2);
}

std::ostream& operator<<(std::ostream& os, const Position& p){
	return os << '(' << p.x << ", " << p.y << ", " << p.z << ')';
}

///\param p the point to which the muon path's closest approach point is to be computed
///\param m a point on the muon track
///\param dir the direction of the muon track
///\pre dir should be a unit vector
Position closestApproach(const Position& p, const Position& m, const Position& dir){
	//We want to solve for point c such that (p-c)*dir = 0
	//let x be the unknown distance along the track from m to c: c=m+x*dir
	//so (p-m-x*dir)*dir = 0
	//p*dir - m*dir - x*dir*dir = 0
	//x*dir*dir = p*dir - m*dir
	double x = (p*dir - m*dir)/dir.mag2();
	return m+x*dir;
}

Position directCherenkovPoint(const Position& p, const Position& m, const Position& dir, double theta_c){
	Position c=closestApproach(p,m,dir);
	double cDist=(p-c).magnitude(); //closest approach distance
	if(cDist==0)
		return c;
	//tan(theta_c) = cDist/x
	//where x is the distance between the point of closest approach and the Cherenkov emission point
	double x = cDist/tan(theta_c);
	Position pc=c-x*dir;
	return pc;
}

int main(int argc, char* argv[]){
	const double c=0.299792458; //m/ns
	const double n=1.34744898746;
	const double ngroup=1.376924748647469;
	const double theta_c=acos(1/n); //rad
	
	double d=100;
	double d_att=2.3*35;
	std::size_t iterations=1000000;
	std::size_t seed=7262;
	unsigned int dimensions=3;
	
	OptionParser op;
	op.addOption({"d","dist","distance"},d,"Distance between modules");
	op.addOption({"d_att"},d_att,"Maximum distance at which photon reception is treated as possible");
	op.addOption({"i","iterations"},iterations,"Number of Monte Carlo trials");
	op.addOption({"s","seed"},seed,"RNG seed");
	op.addOption({"dim"},dimensions,"Number of dimensions to simulate; must be 2 or 3");
	std::vector<std::string> positionals = op.parseArgs(argc, argv);
	if(op.didPrintUsage())
		return 0;
	
	if(dimensions!=2 && dimensions!=3){
		std::cerr << "Number of dimensions must be 2 or 3" << std::endl;
		return 1;
	}
	
	#ifdef USE_PCG
	pcg32 rng(seed);
	#else
	std::mt19937 rng(seed);
	#endif
	
	//let the two modules be at (0,0,0) and (d,0,0)
	Position p1(0,0,0); //first module
	Position p2(d,0,0); //second module
	
	std::uniform_real_distribution<double> posDist(-d_att,d_att);
	double tMin=std::numeric_limits<double>::infinity();
	double tMax=-std::numeric_limits<double>::infinity();
	for(std::size_t i=0; i<iterations; i++){
		//sample two arbitrary points near p1 and p2 for the muon to pass through
		Position disp1=Position(posDist(rng),dimensions==3?posDist(rng):0,posDist(rng));
		Position disp2=Position(posDist(rng),dimensions==3?posDist(rng):0,posDist(rng));
		//Position disp1=Position(50,dimensions==3?posDist(rng):0,posDist(rng));
		//Position disp2=Position(-50,dimensions==3?posDist(rng):0,posDist(rng));
		
		Position mp1=p1+disp1;
		//Position mp1=p1;
		Position mp2=p2+disp2;
		//Position mp2=p2;
		Position muonDir=(mp2-mp1).normed();
		Position pc1=directCherenkovPoint(p1,mp1,muonDir,theta_c);
		Position pc2=directCherenkovPoint(p2,mp1,muonDir,theta_c); //either mp1 or mp2 should be correct here
		
		//skip attempts where the cherenkov emission points ended up beyond the allowed attenuation distance
		if((p1-pc1).mag2()>d_att*d_att){
			i--;
			continue;
		}
		if((p2-pc2).mag2()>d_att*d_att){
			i--;
			continue;
		}
		//return the time at which the muon passes through p, assuming it is on the track
		//take the time at which the muon is at mp1 as the reference point t==0
		auto muonTime=[&](Position p){
			Position v=p-mp1;
			return (v*muonDir)/c;
		};
		
		double t1=muonTime(pc1) //time for muon to travel from pc1 to pc2
		         +(p1-pc1).magnitude()/(c/n); //time for light to travel from pc1 to p1
		double t2=muonTime(pc2) //time for muon to travel from mp1 to pc2
		        +(p2-pc2).magnitude()/(c/n); //time for light to travel from pc2 to p2
		double td=t2-t1;
		if(td>=0 && td<tMin){
			tMin=td;
// 			std::cout << "New tMin: " << tMin << std::endl;
// 			std::cout << "  pc1: " << pc1 << std::endl;
// 			std::cout << "  pc2: " << pc2 << std::endl;
// 			std::cout << "  dc1: " << (p1-pc1).magnitude() << std::endl;
// 			std::cout << "  dc2: " << (p2-pc2).magnitude() << std::endl;
// 			std::cout << "  pc1->pc2: " << (pc1-pc2).magnitude() << std::endl;
// 			std::cout << "  emission angle 1: " << acos(muonDir.unit()*(p1-pc1).unit()) << std::endl;
// 			std::cout << "  emission angle 2: " << acos(muonDir.unit()*(p2-pc2).unit()) << std::endl;
// 			std::cout << "  Light time pc1 to p1: " << (p1-pc1).magnitude()/(c/n) << std::endl;
// 			std::cout << "  Muon time pc1 to pc2: " << (pc2-pc1).magnitude()/c << std::endl;
// 			std::cout << "  Light time pc2 to p2: " << (p2-pc2).magnitude()/(c/n) << std::endl;
		}
		if(td>tMax){
			tMax=td;
// 			std::cout << "New tMax: " << tMax << std::endl;
// 			std::cout << "  pc1: " << pc1 << std::endl;
// 			std::cout << "  pc2: " << pc2 << std::endl;
// 			std::cout << "  dc1: " << (p1-pc1).magnitude() << std::endl;
// 			std::cout << "  dc2: " << (p2-pc2).magnitude() << std::endl;
// 			std::cout << "  pc1->pc2: " << (pc1-pc2).magnitude() << std::endl;
// 			std::cout << "  muon dir: " << muonDir << std::endl;
// 			std::cout << "  emission angle 1: " << acos(muonDir.unit()*(p1-pc1).unit()) << std::endl;
// 			std::cout << "  emission angle 2: " << acos(muonDir.unit()*(p2-pc2).unit()) << std::endl;
// 			double alpha=acos((p1-p2)*(pc2-p2)/((p1-p2).magnitude()*(pc2-p2).magnitude()));
// 			std::cout << "  alpha: " << alpha << std::endl;
// 			std::cout << "  cos(alpha): " << cos(alpha) << std::endl;
// 			std::cout << "  n*cos(alpha+theta_c): " << n*cos(3.14159265-alpha-theta_c) << std::endl;
// 			std::cout << "  Light time pc1 to p1: " << (p1-pc1).magnitude()/(c/n) << std::endl;
// 			std::cout << "  Muon time pc1 to pc2: " << (pc2-pc1).magnitude()/c << std::endl;
// 			std::cout << "  Light time pc2 to p2: " << (p2-pc2).magnitude()/(c/n) << std::endl;
		}
	}
	std::cout << "Observed time difference range: [" <<  << ", " << tMax << ']' << std::endl;
	
	double sc=sin(theta_c);
	
// 	std::cout << "Min intermediate: " << sqrt(d*d-4*d_att*d_att*sc*sc)/c << std::endl;
// 	
// 	//this is valid as the minimum possible time difference for distnces d>2*d_att
// 	double t_min_ana=(d_att/c)*(sqrt(1-sc*sc)+sqrt((d*d)/(d_att*d_att)-sc*sc)-n);
// 	double t_max_ana=(d_att/c)*(-sqrt(1-sc*sc)+sqrt((d*d)/(d_att*d_att)-sc*sc)+n);
// 	std::cout << "Analytic result: [" << t_min_ana << ", " << t_max_ana << ']' << std::endl;
	
// 	double D=(1/n)-n;
// 	double cross=(-d_att/(2*D))*sqrt(9*pow(sc,4)+10*sc*sc*D*D+pow(D,4));
// 	std::cout << "Cross-over distance: " << cross << std::endl;
 	
//  	double alpha=0;//1.14665;
//  	std::cout << "max close:  " << (d/(c*sin(theta_c)))*(n*sin(alpha+theta_c)-sin(alpha)) << std::endl;
//   	std::cout << "max close: " << n*d/c << std::endl;

}