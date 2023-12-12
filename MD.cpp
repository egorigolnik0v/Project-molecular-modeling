#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

template <typename type>
class Particle {

public:

type x;
type y;
type z;
type vx;
type vy;
type vz;
type ax;
type ay;
type az;

Particle(){}


Particle(type x, type y, type z, type vx, type vy, type vz, type ax, type ay, type az) : x(x),y(y),z(z),vx(vx),vy(vy),vz(vz),ax(ax),ay(ay),az(az){}

void move(float dt, float limit){

vx  += ax * dt;
vy += ay * dt;
vz += az * dt;
x += vx * dt;
y += vy * dt;
z += vz * dt;

if (x  > limit) {x = -limit; }
if (y  > limit) {y = -limit; }
if (z  > limit) {z = -limit; }
if (x  < -limit) {x = limit; }
if (y  < -limit) {y = limit; }
if (z  < -limit) {z = limit; }
}

};




int main() {

	int k;
	cout << "Molecular dynamic" << endl;
	cout << "Potential: Lennard-Jones" << endl;
	cout << "Number of particles per side of the cube : ";
	cin >> k;
	int N = k * k * k;
	cout << "Number of particles: " << N << endl;
	long double limit = k / 2 + 0.5;
	cout << "recommended: " << limit*2 << " (σ)   length of simulation cell:  (σ) ";
	cin >> limit;
	limit = limit/2;
	cout << "recommended: " << limit << " (σ)  cut LJ-potential after:  (σ) ";
	long double radius = limit;
	cin >> radius;
	long double radius_2 = radius * radius;
	int step = 51;
	long double* V = new long double[N];
	long double* list = new long double[step];
	long double* ax_ = new long double[N];
	long double* ay_ = new long double[N];
	long double* az_ = new long double[N];
	int t = 0;
	float interval = 1;
	cout << "recommended: " << interval << "  (σ)   distance between particles:  (σ) ";
	cin >> interval;
	float dt = 0.001;
	cout << "recommended: " << dt << "   time of one iteration step: ";
	cin >> dt;
	cout << "Number of iteration steps: ";
	int iteration = 1;
	cin >> iteration;
	long double V_max;
	long double Potential_Energy;
	long double Kenetic_Energy;
	string line;

	cout << "Setting start parametres: ";
	
	vector<Particle<float>> Particles;
	
	for (int i = 0; i < N; i++) {

		Particle<float> A = Particle<float>();

		A.x = interval * ((i % k) - int(k / 2)) + 0.5;
		A.y = interval * (((i / k) % k) - int(k / 2) + 0.5);
		A.z = interval * ((i / (k * k)) - int(k / 2)) + 0.5;
		int sign_x = rand();
		float v_x = rand();
		v_x = v_x / RAND_MAX;
		sign_x = (sign_x % 3) - 1;
		A.vx = v_x * sign_x;

		int sign_y = rand();
		float v_y = rand();
		v_y = v_y / RAND_MAX;
		sign_y = (sign_y % 3) - 1;
		A.vy = v_y * sign_y;

		int sign_z = rand();
		float v_z = rand();
		v_z = v_z / RAND_MAX;
		sign_z = (sign_z % 3) - 1;
		A.vz = v_z * sign_z;

		A.az = 0;
		A.ay = 0;
		A.az = 0;
		
		Particles.push_back(A);
	}
	cout << "Done" << endl;

	for (int i = 0; i < step; i++) { list[i] = 0; }

	ofstream f("Data.txt");
	ofstream e("Energy.txt");
	ofstream K("Energy_K.txt");
	ofstream l("List.txt");
	ofstream p("Data_forces.txt");
	ifstream in("positions_lammps");
	ifstream in_2("forces_lammps");

	bool choice;
	cout << "Set coordinates using positions_lammps ? [yes/no]  : ";
	string A;
	cin >> A;
	if (A=="yes") {choice = true;}
	else{choice = false;}

	int q = 0;
	if (in.is_open() and choice) {

		int m = 0;
		while (getline(in, line)) {
			int t = 1;
			q++;
			if (q>9){
			stringstream ss(line);
			string word;
			while (ss >> word and q>9) {
				if (t == 1) { Particles[m].x = stof(word); }
				if (t == 2) { Particles[m].y = stof(word); }
				if (t == 3) { Particles[m].z = stof(word); }
				t++;
			}m++;}
		}
	}
	
	bool choice2;
	cout << "Set forces using forces_lammps ? [yes/no]  : ";
	string B;
	cin >> B;
	if (B=="yes") {choice2 = true;}
	else{choice2 = false;}
	
	int r = 0;
	if (in_2.is_open() and choice2) {

		int m = 0;
		while (getline(in_2, line)) {
			int t = 1;
			r++;
			if (r>9){
			stringstream ss(line);
			string word;
			while (ss >> word) {
				if (t == 1) { ax_[m] = stof(word); }
				if (t == 2) { ay_[m] = stof(word); }
				if (t == 3) { az_[m] = stof(word); }
				t++;
			}m++;}
		}
	}
	
	cout << "Start modeling ...";

	for (int i1 = 0; i1 < iteration; i1++) {

		Potential_Energy = 0;
		Kenetic_Energy = 0;
		V_max = 0;

		for (int i = 0; i < N; i++) {

			Particles[i].ax = 0;
			Particles[i].ay = 0;
			Particles[i].az = 0;
		}

		for (int i = 0; i < (N - 1); i++) {
			for (int j = i; j < N; j++) {
				if (i != j) {
					long double dist_ij_x = Particles[i].x - Particles[j].x;
					long double dist_ij_y = Particles[i].y - Particles[j].y;
					long double dist_ij_z = Particles[i].z - Particles[j].z;

					if (dist_ij_x > limit) { dist_ij_x = dist_ij_x - 2 * limit; }
					if (dist_ij_y > limit) { dist_ij_y = dist_ij_y - 2 * limit; }
					if (dist_ij_z > limit) { dist_ij_z = dist_ij_z - 2 * limit; }
					if (dist_ij_x < -limit) { dist_ij_x = dist_ij_x + 2 * limit; }
					if (dist_ij_y < -limit) { dist_ij_y = dist_ij_y + 2 * limit; }
					if (dist_ij_z < -limit) { dist_ij_z = dist_ij_z + 2 * limit; }

					long double dist_ij_2 = dist_ij_x * dist_ij_x + dist_ij_y * dist_ij_y + dist_ij_z * dist_ij_z;
					long double dist_ij__2 = 1 / dist_ij_2;
					long double dist_ij__6 = dist_ij__2 * dist_ij__2 * dist_ij__2;
					long double dist_ij__12 = dist_ij__6 * dist_ij__6;
					long double dist_ij__8 = dist_ij__6 * dist_ij__2;
					long double dist_ij__14 = dist_ij__12 * dist_ij__2;
					long double forc_ij__r = 4 * (12 * dist_ij__14 - 6 * dist_ij__8);
					if (dist_ij_2 <= radius_2) {

						Particles[i].ax += dist_ij_x * forc_ij__r;
						Particles[i].ay += dist_ij_y * forc_ij__r;
						Particles[i].az += dist_ij_z * forc_ij__r;
						Particles[j].ax -= dist_ij_x * forc_ij__r;
						Particles[j].ay -= dist_ij_y * forc_ij__r;
						Particles[j].az -= dist_ij_z * forc_ij__r;

						Potential_Energy += 4 * (1 * dist_ij__12 - 1 * dist_ij__6);
					}
				}
			}
		}

		ofstream("Data_forces.txt", ios::app);
		for (int i = 0; i < N; i++) {
			p << setprecision(15) << abs(Particles[i].ax / ax_[i]) << "\t" << abs(Particles[i].ay / ay_[i]) << "\t" << abs(Particles[i].az / az_[i]) << endl;
		}

		for (int i = 0; i < N; i++) {

			Particles[i].move(dt, limit);

			V[i] = Particles[i].vx * Particles[i].vx + Particles[i].vy * Particles[i].vy + Particles[i].vz * Particles[i].vz;
			Kenetic_Energy = Kenetic_Energy +
				0.5 * V[i];

			if (V[i] > V_max) { V_max = V[i]; }
			
			
		}

		for (int j = 0; j < step; j++) {
			for (int i = 0; i < N; i++) {
				if (abs(V[i] - V_max * j / 50) <= V_max / 100) { list[j]++; }
			}
		}

		if (0 == 0) {

			ofstream("List.txt", ios::app);
			for (int i = 0; i < step; i++) {
				if (i < (step - 1)) { l << list[i] << ' '; }
				else { l << list[i] << " " << V_max << " " << Kenetic_Energy << endl; }
				list[i] = 0;
			}

			ofstream("Data.txt", ios::app);
			f << N << endl;
			f << endl;
			for (int i = 0; i < N; i++) {
				f << Particles[i].x << " " <<Particles[i].y << " " <<Particles[i].z << endl;
			}

			ofstream("Energy.txt", ios::app);

			e << (Kenetic_Energy + Potential_Energy) << endl;

			ofstream("Energy_K.txt", ios::app);
			K << Kenetic_Energy << endl;
		}

		t = t + 1;
	}

	f.close();
	e.close();
	K.close();
	l.close();
	in.close();
	cout << " Complited" << endl;
	return 0;
}
