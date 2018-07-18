#include "duc.h"

#if 1
//______________________________________________________________________________

void duc(const std::array<DATA_T, L_INPUT>&  din_i, const std::array<DATA_T, L_INPUT>&  din_q,
        DATA_T dout[L_OUTPUT], const incr_t incr) {

    //#pragma HLS interface ap_ctrl_none port=return  // to avoid bubble
    //#pragma HLS interface s_axilite port=incr
    //#pragma HLS interface axis depth=300  port=din_i
    //#pragma HLS interface axis depth=300  port=din_q
    //#pragma HLS interface axis depth=4800 port=dout
    #pragma HLS interface axis depth=L_INPUT  port=din_i
    #pragma HLS interface axis depth=L_INPUT  port=din_q
    #pragma HLS interface axis depth=L_OUTPUT port=dout
    
    #pragma HLS interface ap_stable port=incr
    
    #pragma HLS dataflow
    
    static filterStageClassTwoChan <L_INPUT> f0;
    static dds_class <L_OUTPUT, dds_t, DATA_T> dds_0;
    
    std::array<DATA_T, L_OUTPUT>  dout_i;
    std::array<DATA_T, L_OUTPUT>  dout_q;
    dds_t  dds_cos[L_OUTPUT];
    dds_t  dds_sin[L_OUTPUT];
    //DATA_T  dout_ii[L_OUTPUT];
    //DATA_T  dout_qq[L_OUTPUT];

    f0.process(din_i, dout_i, din_q, dout_q);
    /*for (auto i = 0; i < dout_q.size(); ++i) {
        dout_qq[i] = dout_q[i];
        dout_ii[i] = dout_i[i];
    }*/
    dds_0.process_frame(incr, dds_cos, dds_sin);
    dds_0.mix(dds_cos, dds_sin, dout_i.data(), dout_q.data(), dout);

}

/*void filters(const std::array<DATA_T, L_INPUT>&  din_i, std::array<DATA_T, L_OUTPUT>&  dout_i,
             const std::array<DATA_T, L_INPUT>&  din_q, std::array<DATA_T, L_OUTPUT>&  dout_q) {

    #pragma HLS interface axis depth=L_INPUT  port=din_i
    #pragma HLS interface axis depth=L_INPUT  port=din_q
    #pragma HLS interface axis depth=L_OUTPUT port=dout_i
    #pragma HLS interface axis depth=L_OUTPUT port=dout_q
    #pragma HLS dataflow
    
    static filterStageClassTwoChan <L_INPUT> f0;
    f0.process(din_i, dout_i, din_q, dout_q);
}


void process_single_filter(const std::array<DATA_T, L_INPUT>& din, std::array<DATA_T, 2*L_INPUT>& dout) {
#pragma HLS INLINE
#pragma HLS dataflow
    static filterStageClassTwoChan <L_INPUT> f0;
    f0.process_single_fir(din, dout);
}
*/
#endif





