#include <ctime> /* clock_t, clock, CLOCKS_PER_SEC */

class Clock {
	
	public:
	
	// Constructor
	Clock();
	
	// Start and stop the time clock
	void start() {t_start = clock();};
	void stop() {t_stop = clock();};
	
	// Print the time
	void print();
	
	private:
	
	// Time clock
	clock_t t_start, t_stop;
		
};