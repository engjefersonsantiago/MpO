#include "duc.h"
#include <math.h>
#include <iterator>
#include <array>
#include <algorithm>


//______________________________________________________________________________
void duc_sw ( DATA_T din_i[L_INPUT], 
	          DATA_T din_q[L_INPUT], 
		      DATA_T dout[L_OUTPUT],
              incr_t incr ) {

static filterStageClassTwoChan <L_INPUT> f0;
static dds_class <L_OUTPUT, dds_t, DATA_T> dds_0;

DATA_T dout_i[L_OUTPUT];
DATA_T dout_q[L_OUTPUT];
dds_t  dds_cos[L_OUTPUT];
dds_t  dds_sin[L_OUTPUT];

std::array<DATA_T, L_INPUT> din_ii;
std::array<DATA_T, L_OUTPUT> dout_ii;
std::array<DATA_T, L_INPUT> din_qq;
std::array<DATA_T, L_OUTPUT> dout_qq;

for (auto i = 0; i < din_ii.size(); ++i) {
    din_ii[i] = din_i[i];
    din_qq[i] = din_q[i];
}
f0.process(din_ii, dout_ii, din_qq, dout_qq);
dds_0.process_frame(incr, dds_cos, dds_sin);
     for (auto i = 0; i < dout_qq.size(); ++i) {
        dout_q[i] = dout_qq[i];
        dout_i[i] = dout_ii[i];
    }

dds_0.mix(dds_cos, dds_sin, dout_i, dout_q, dout);

}


//______________________________________________________________________________
int main() {

  ofstream fp_result("out.dat",ofstream::out); 
  fp_result.precision(25); 
  fp_result.width(30);

    double tmp;

    const int nofset = 10;

    DATA_T din_i[nofset][L_INPUT];
    DATA_T din_q[nofset][L_INPUT];
    DATA_T dout[nofset][L_OUTPUT];
    DATA_T swout[nofset][L_OUTPUT];

incr_t  incr;
//incr = 303.33333*pow(2,31-11);
incr = 303.33333*pow(2,-11);

    int jj = 0;
	float diff;
	int err_cnt = 0;

    for (int i=0; i<nofset; i++) {
        
        for (int k=0; k<L_INPUT; k++) {
            tmp = cos(2*M_PI*(0.5+(double)jj*33)/(1024));
            din_i[i][k] = tmp;

#if 0
            tmp = 0.3*cos(2*M_PI*(0.5+(double)jj*10)/(1024));
            din_q[i][k] = tmp;
#endif
            tmp = sin(2*M_PI*(0.5+(double)jj*33)/(1024));
            din_q[i][k] = tmp;

            jj++;
        }
        std::array<DATA_T, L_INPUT> din_ii;
        std::array<DATA_T, L_OUTPUT> dout_ii;
        std::array<DATA_T, L_INPUT> din_qq;
        std::array<DATA_T, L_OUTPUT> dout_qq;
        
        for (auto ii = 0; ii < din_ii.size(); ++ii) {
            din_ii[ii] = din_i[i][ii];
            din_qq[ii] = din_q[i][ii];
        }
        duc(din_ii, din_qq, dout[i], incr);
        duc_sw(din_i[i], din_q[i], swout[i], incr);

        //for (int k=0; k<L_OUTPUT; k++) 
            //fp_result << dout[i][k] << endl;

        // compare with expected
        for (int k=0; k<L_OUTPUT; k++)  {
		  diff = dout[i][k] - swout[i][k];
          fp_result << dout[i][k] <<", "<<swout[i][k] << endl;
		  if (diff != 0) err_cnt++;
	    }

    }

  if (diff==0) {
      cout <<"============== Test passed "<<endl;
      return 0;
  } else {
      cout <<"============== Test failed, error count =  "<< err_cnt << endl;
      return 1;
  }


}



