#include <iostream>
using namespace std;

#include "clock.h"

Clock::Clock() {
	t_start = 0;
	t_stop  = 0;
}

void Clock::print() {

	if (t_start == 0 || t_stop == 0) {
		cout << "WARNING: please, start and stop the Clock before print." << endl;
		return;
	}
	
	cout << endl;
	cout << "Starting clock: " << t_start << endl;
	cout << "Ending clock: " << t_stop << endl;
	double elapsed_time = t_stop - t_start;	
	cout << "Time for this run: " << elapsed_time << " clocks =  ";
	cout << (double)elapsed_time/CLOCKS_PER_SEC << " seconds";
	cout << endl << endl;
}