#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allocore/io/al_App.hpp"
using namespace al;

struct Graph {
  
  Graph() {}
  
  void draw(Graphics& g, Mesh& m) {
  m.primitive(Graphics::LINES);
  m.vertex(0,0,0);
  m.vertex(3,0,0);
  m.vertex(0,0.05,0);
  m.vertex(0,-0.05,0);
  m.vertex(1,0.05,0);
  m.vertex(1,-0.05,0);
  m.vertex(2,0.05,0);
  m.vertex(2,-0.05,0);
  m.vertex(3,0.05,0);
  m.vertex(3,-0.05,0);

  for (unsigned i = 0; i < 36; i++) {
    m.vertex(i/12.0f,0.01,0);
    m.vertex(i/12.0f,-0.01,0);
  }
  g.draw(m);
  }
};

struct AlloApp : App {

int iBin;
int ind=0;
int numpartials;
int i;
float ampl[1024];
float frequ[1024];
float diss[400];
Graph axis;
Mesh graph;

  AlloApp() {
    nav().pos(0,0,10);
    initWindow();

	  printf("\n");
	  printf("Please enter the number of partials and press enter:\n");
	  scanf( "%d", &numpartials );
	  printf("\n");
	
	

	  for(i=1; i<=numpartials; i++) {
		  printf("Please enter the frequency of partial %d and press enter:\n", i);
		  printf("\n");
		  scanf("%f", &frequ[i]);
		  printf("\n");
		  printf("Please enter the amplitude of partial %d and press enter:\n", i);
		  printf("\n");
		  scanf("%f", &ampl[i]);
	  }
	
	  //dissonance calculation part is as follows;
	
		float dstar= (float) (0.24);  // this is the point of maximum dissonance - the value is derived from a model
							// of the Plomp-Levelt dissonance  curves for all frequencies.
		float s;
		float s1 = (float) (0.0207) ;// s1 and s2 are used to allow a single functional form to interpolate beween
		float s2 = (float) (18.96);  //  the various P&L curves of different frequencies by sliding, stretching/compressing
		       				//  the curve so that its max dissonance occurse at dstar. A least-square-fit was made
							// to determine the values.
		
		float c1 = 5;       // these parameters have values to fit the experimental data of Plomp and Levelt
		float c2 = -5;
		float a1 = (float) (-3.51); // theses values determine the rates at which the function rises and falls and 
		float a2 = (float) (-5.75); // and are based on a gradient minimisation of the squared error between 
						  //  Plomp and Levelt's averaged data and the curve
		
		// If the point of maximum dissonance for a base frequency f occurs at dstar, then the dissonance between
		// f1 with amp1 and f2 with amp2, is - amp1*amp2(e^-a1s[f2-f1] - e^-a2s[f2-f1]) where s = dstar/(s1f1+f2).....
		
		float fdif; // will hold the difference between frequency values


		float arg1;
		float arg2;
		float dnew;
		float d =0;
		float exp1;
		float exp2;
		
		ind = 0;
		int i;
		int j;
		int k;
				
		float fmin; //lowest partial - i.e, f1  

		float lowint = 1 ;
		float highint = 4 ; // sets the upper limit of intervals for which dissonance will be calculated.
		float interval;
		float inc = (float) (0.01);
		
		float size=((highint-lowint)/inc);
		
	//	float diss[400]; // declare array of dissonance values for each step from lowint to highint.(1-2)
		float intervals[400]; //declare array of intervals for final curve
		float allpartialsatinterval[1024]; // declare array of partials pitched up by 'interval'
		
		for (interval = lowint; interval <= highint; interval += inc) 
		//perform below for interval values lowint - highint.
			{
			d=0;
				for (k = 1; k <= numpartials; k++) 
					// fill allpartialsatinterval array with each element of freq array multiplied by interval.
					
					{
						allpartialsatinterval[k] = interval*(frequ[k]);// what is in these arrays first time round?
					}
				// Calculate the dissonance between frequ[] and interval*frequ[]
				for (i=1; i<=numpartials; i++)
					{
						for (j=1; j<=numpartials; j++)
						{
							
								// if an element of allpartials is less than this element of freq
								// then give its value to fmin
							if (allpartialsatinterval[j]<frequ[i])//if 1.x by each fq component is less than the current component
								{

									fmin = allpartialsatinterval[j] ;//fmin takes on the lesser of freq*1.x and frequ
								}
							else 
								// fmin takes on the current frequ value
								{  
									fmin = frequ[i] ;
								}
							
							s=dstar/(s1*fmin+s2); // define s with interpolating values s1 and s2.
									
							fdif=(float) (fabs(allpartialsatinterval[j]-frequ[i])) ; // (fabs gives absolute value)
																		//	establishes the frequency difference					
								
							
							arg1=a1*s*fdif; arg2=a2*s*fdif ; // gives arg1/arg2 the powers for e [i.e., as(f2-f1)] 
															// in the equation for dissonance given above.
							
						
							    exp1=(float) (exp(arg1)); // EXP returns a value containing e
														//(the base of natural logarithms) raised to the specified power.
						 
							
								exp2=(float) (exp(arg2));
						
								
							if (ampl[i]<ampl[j])
								{

								dnew=ampl[i]*(c1*exp1+c2*exp2);//  The lesser of amp(i) and amp(j) is used
																//the idea is to give larger dissonance for louder sounds 
								}								// and to smoothly go to zero as one of the partials disappears 
							else{
								dnew=ampl[j]*(c1*exp1+c2*exp2); // use ampl(j) if it is smaller
								
								}
							d = d+dnew;  // keep adding the dissonances of each loop, where d iterates the dissonance of 1 partial(frequ[i]) with 
										// each(all) of the timbre's other partials as they'd be at the current interval. This is done for 
										// for each partial(freq[i]) and every interval.
						}				
					}		
						
					diss[ind] = d;//fill  diss array with dissonances
					intervals[ind] = interval; // fill  intervals array with intervals 
					printf("ind is: %d .\n", ind);
					printf("Interval %f - Dissonance %f \n",intervals[ind], diss[ind]);					
					ind++;
				}
			
		int dissvalues = (int) ((highint-lowint)/inc);	
		printf("There are %d dissonance values\n", dissvalues);	
  }
		
  void onDraw(Graphics& g) {
    Mesh m;
    m.primitive(Graphics::LINE_STRIP);

    for (int i = 0; i < 300; i++)
      m.vertex(i/100.0f, diss[i], 0);
      
    axis.draw(g,graph);

    g.draw(m);
  }
};

int main() {AlloApp().start();}
