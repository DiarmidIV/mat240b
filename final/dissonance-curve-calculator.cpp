#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "allocore/io/al_App.hpp"
#include "allocore/math/al_Ray.hpp"

#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Envelope.h"
#include "Gamma/Gen.h"

#include "alloGLV/al_ControlGLV.hpp"
#include "GLV/glv.h"

using namespace al;


struct CalculateCurve {
  int ind;
  int numpartials = 8;
  int i;
  float ampl[1024];
  float frequ[1024];
  float diss[400];
  float intervals[400];

  float dstar = (float) (0.24);  
  float s;
  float s1 = (float) (0.0207);
  float s2 = (float) (18.96);
  float c1 = 5;
  float c2 = -5;
  float a1 = (float) (-3.51); 
  float a2 = (float) (-5.75); 
  float fdif; 
  float arg1;
  float arg2;
  float dnew;
  float d = 0;
  float exp1;
  float exp2;
  int j;
  int k;
  float fmin;  
  float lowint = 1;
  float highint = 4;
  float interval;
  float inc = (float) (0.01);
  float size=((highint-lowint)/inc);
  float allpartialsatinterval[1024];  

  CalculateCurve() {}

  void calculate(float freqs[], float amps[]) {
    ind = 0;
    for (int m = 0; m < numpartials; m++) {
    frequ[m] = freqs[m];
    ampl[m] = amps[m];
    }  
		for (interval = lowint; interval <= highint; interval += inc)	{
      d = 0;
			for (k = 1; k <= numpartials; k++) {
        allpartialsatinterval[k] = interval*(frequ[k]);
			}
			for (i = 1; i <= numpartials; i++) {
				for (j = 1; j <= numpartials; j++) {
					if (allpartialsatinterval[j] < frequ[i]) {
						fmin = allpartialsatinterval[j];
								}
						else {  
									fmin = frequ[i];
								}
							
            s = dstar / (s1 * fmin + s2); 			
						fdif = (float) (fabs(allpartialsatinterval[j] - frequ[i])); 
						arg1 = a1 * s * fdif; arg2 = a2 * s * fdif;  
						exp1 = (float) (exp(arg1));
						exp2 = (float) (exp(arg2));
								
						if (ampl[i] < ampl[j]) {
						  dnew = ampl[i] * (c1 * exp1 + c2 * exp2);
            }								 
            else {
              dnew = ampl[j] * (c1 * exp1 + c2 * exp2);
            }
            d = d + dnew; 
        }				
      }		
						
      diss[ind] = d;
			intervals[ind] = interval;  
		//	printf("ind is: %d .\n", ind);
		//	printf("Interval %f - Dissonance %f \n",intervals[ind], diss[ind]);					
			ind++;
    }
			
		int dissvalues = (int) ((highint-lowint)/inc);	
	//	printf("There are %d dissonance values\n", dissvalues);	
  }
};

struct Data {
  float interval;
  float dissonance;

  Data(float initInterval, float initDissonance) {
    interval = initInterval;
    dissonance = initDissonance;
  }
};

bool compareInterval (Data a, Data b) { return (a.interval < b.interval); } 

bool compareDissonance (Data a, Data b) { return (a.dissonance < b.dissonance); }

struct Graph {
  Mesh m;

  Graph() {
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
  }
  
  void draw(Graphics& g) {
    g.color(1,1,1,1);
    g.draw(m);
  }
};

struct ScaleDegree {
  Vec3f position;
  float transposition;
  bool state = false;
  Color active = {1, 0, 0, 1};
  Color inactive = {0, 0, 1, 1,};
  Color currentColor = {inactive};
  gam::Sine<> sine[8];
  gam::Seg<> env;

  ScaleDegree(Vec3f initPosition) {
    position = initPosition;
    transposition = position.x + 1.0f;
    env = {0,0,0,0};
      
    env.length(0.05);
  }

  void colorState() {
    if (state == false) {
      currentColor = active;
      state = true;
    }
    else if (state == true) {
      currentColor = inactive;
      state = false;
    }
  }
  
  float audio(float freqs[8], float amps[8]) {
    float s = 0;
    for (int i = 0; i < 8; i++) {
      sine[i].freq(freqs[i] * transposition);
      s += sine[i]() * amps[i];
    }
    env = state;
    s *= env();
    return s;
  }

  void draw(Graphics& g, Mesh& m) {
    g.pushMatrix();
    g.translate(position);
    g.color(currentColor);
    g.draw(m);
    g.popMatrix();
  }
};

struct AlloApp : App {

CalculateCurve calc;  
std::vector <Data> data;
std::vector <ScaleDegree> scaleDegree;
Graph axis;
Mesh graph;
Mesh sphere;
float radius = 0.01f;

gam::Seg<> env;

GLVBinding gui;
glv::Slider slider;
glv::Sliders sliders;
glv::Sliders sliders2;
float f[8] = {0,0,0,0,0,0,0,0};
float a[8] = {0,0,0,0,0,0,0,0};
bool redrawGraph = false;

glv::Table layout;

  AlloApp() {
    nav().pos(0,0,10);
    navControl().useMouse(false);
   
    initWindow();
    
    addSphere(sphere, radius);
      
    env.length(0.05);

    calc.calculate(f,a);

    for (unsigned i = 1; i < 299; i++) {
      if ( calc.diss[i] < calc.diss[i-1])
        if (calc.diss[i] < calc.diss[i+1]) {
          Data d = {calc.intervals[i],calc.diss[i]}; 
          data.push_back(d); 
        }
    }
    
    std::sort (data.begin(), data.end(), compareDissonance);
    Vec3f position = {0,0,0};
    ScaleDegree d = {position};
    scaleDegree.push_back(d);
    for (unsigned i = 0; i < data.size(); i++) {
      Vec3f position = {(data[i].interval-1.0f),0,0};
      ScaleDegree d = {position};
      scaleDegree.push_back(d);
      std::cout << data[i].interval << ' ' << data[i].dissonance << std::endl;
    }

    scaleDegree[0].colorState();

    gui.bindTo(window());
    gui.style().color.set(glv::Color(0.7), 0.5);
    layout.arrangement("x x");
    layout << new glv::Label("frequency");
    layout << new glv::Label("amplitude");
    sliders = { glv::Rect(100, 64), 1, 8, 1 }; 
    layout << sliders;
    sliders2 = { glv::Rect(100, 64), 1, 8, 1 }; 
    layout << sliders2;
    layout.arrange();
    gui << layout;
  
    initAudio();
  }
	
  void onAnimate(double dt) {
    
    for (int i = 0; i < 8; i++) {
      f[i] = sliders.getValue(i) * 1000;
      a[i] = sliders2.getValue(i);
    }
    
    if (redrawGraph == true) {
      calc.calculate(f,a);
      data.clear();   
      scaleDegree.clear();

      for (int i = 0; i < 8; i++) {
        std::cout << "FREQUENCY: " << f[i] << " " << "AMPLITUDE: " << a[i] << std::endl;
      }

      for (unsigned i = 1; i < 299; i++) {
        if (calc.diss[i] < calc.diss[i-1])
          if (calc.diss[i] < calc.diss[i+1]) {
            Data d = {calc.intervals[i],calc.diss[i]}; 
            data.push_back(d); 
          }
        }
      
      std::cout << "MINIMA: " << std::endl;
      std::sort (data.begin(), data.end(), compareDissonance);
      Vec3f position = {0,0,0};
      ScaleDegree d = {position};         
      scaleDegree.push_back(d);
      for (unsigned i = 0; i < data.size(); i++) {
        Vec3f position = {(data[i].interval-1.0f),0,0};
        ScaleDegree d = {position}; 
        scaleDegree.push_back(d);
        std::cout << "INTERVAL: " << data[i].interval << ' ' << "DISSONANCE: " << data[i].dissonance << std::endl;   
      }
      
      scaleDegree[0].colorState();

      redrawGraph = false;
    }
  }

  void onDraw(Graphics& g) {
    Mesh m;
    m.primitive(Graphics::LINE_STRIP);

    for (int i = 0; i < 300; i++)
      m.vertex(i/100.0f, calc.diss[i], 0);
      
    for (int i = 0; i < scaleDegree.size(); i++) {
      scaleDegree[i].draw(g,sphere);
    }

    axis.draw(g);

    g.draw(m);
  }

  virtual void onSound(al::AudioIOData& io) {
    gam::Sync::master().spu(audioIO().fps());
    while (io()) {
     
      float s = 0;
      float attenuate = 0;  
      for (int i = 0; i < scaleDegree.size(); i ++) {
        s += scaleDegree[i].audio(f, a);    
        attenuate += scaleDegree[i].state * 8;
       } 

       if (attenuate <= 0) attenuate = 1;
       env = attenuate;
       s /= env(); 
       io.out(0) = s;
       io.out(1) = s;
    
    }
  }

  void onKeyDown(const ViewpointWindow&, const Keyboard& k) {  
    redrawGraph = true;
  }

  virtual void onMouseDown(const ViewpointWindow& w, const Mouse& mouse) {
    Rayd mouse_ray = getPickRay(w, mouse.x(), mouse.y());
    for (int i = 0; i < scaleDegree.size(); i ++) {    
      if (mouse_ray.intersectsSphere(scaleDegree[i].position,radius*1.5)) {
        scaleDegree[i].colorState();
      }
    }
  }

};

int main() {AlloApp().start();}
