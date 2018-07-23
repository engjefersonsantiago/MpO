#include "duc.h"

#if 1
//______________________________________________________________________________

void duc(const std::array<DATA_T, L_INPUT>&  din_i, const std::array<DATA_T, L_INPUT>&  din_q,
        DATA_T dout[L_OUTPUT], const incr_t incr) {

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

    f0.process(din_i, dout_i, din_q, dout_q);
    dds_0.process_frame(incr, dds_cos, dds_sin);
    dds_0.mix(dds_cos, dds_sin, dout_i.data(), dout_q.data(), dout);

}

#endif





