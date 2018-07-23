#ifndef DUC_H
#define DUC_H

#include <fstream>
#include <iostream>
#include <array>

static constexpr auto L_INPUT  = 200;
static constexpr auto L_OUTPUT = 16*L_INPUT;

#include "dds.h"
#include "fir.h"


//______________________________________________________________________________
// filter spec
// - putting these inside the duc class makes it not compile

// SRRC
static constexpr auto Lsrrc_WHOLE = 65;
const COEF_T cin_srrc[Lsrrc_WHOLE] = {
      #include "srrc_0db.inc"
};

// HB1
static constexpr auto Lhb1_WHOLE   = 23;
const COEF_T cin_hb1[Lhb1_WHOLE] = {
      #include "b1_fp2.inc"
};

// HB2
static constexpr auto Lhb2_WHOLE   = 11;
const COEF_T cin_hb2[Lhb2_WHOLE] = {
      #include "b2_fp2.inc"
};

// HB3
static constexpr auto Lhb3_WHOLE   = 7;
const COEF_T cin_hb3[Lhb3_WHOLE] = {
      #include "b3_fp2.inc"
};
  
static constexpr auto II_SRRC = 16;
static constexpr auto II_HB1  = 8;
static constexpr auto II_HB2  = 4;
static constexpr auto II_HB3  = 2;

//______________________________________________________________________________
void duc(const std::array<DATA_T, L_INPUT>&  din_i, const std::array<DATA_T, L_INPUT>&  din_q,
        DATA_T dout[L_OUTPUT], const incr_t incr);


//______________________________________________________________________________
// multi stage filter
// - 4 stages, SRRC -> HB1-> HB2-> HB3 -> to DAC
//______________________________________________________________________________
template<int l_INPUT> class filterStageClass {

    public:

        FirInterp2<Lsrrc_WHOLE,   l_INPUT,   II_SRRC> srrc;
        FirInterp2Hb<Lhb1_WHOLE, 2*l_INPUT, II_HB1>  hb1;
        FirInterp2Hb<Lhb2_WHOLE, 4*l_INPUT, II_HB2>  hb2;
        FirInterp2Hb<Lhb3_WHOLE, 8*l_INPUT, II_HB3>  hb3;


        //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
        void process(const std::array<DATA_T, l_INPUT>& din, std::array<DATA_T, 16*l_INPUT>& dout) {
        #pragma HLS INLINE
        #pragma HLS dataflow
        
            std::array<DATA_T, 2*l_INPUT> srrc_dout;
            std::array<DATA_T, 4*l_INPUT> hb1_dout;
            std::array<DATA_T, 8*l_INPUT> hb2_dout;
                
            srrc.process_frame(din, srrc_dout);
            hb1.process_frame(srrc_dout, hb1_dout);
            hb2.process_frame(hb1_dout, hb2_dout);
            hb3.process_frame(hb2_dout, dout);
        }


    //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
    //constructor
    filterStageClass() :  
            srrc(cin_srrc),
            hb1(cin_hb1),
            hb2(cin_hb2),
            hb3(cin_hb3)
    {}

};


//______________________________________________________________________________
// 2 channel multi-stage filter
//______________________________________________________________________________
//template<int l_INPUT>
template<int l_INPUT, int l_OUTPUT = 16*l_INPUT> class filterStageClassTwoChan {

    private:
        filterStageClass<l_INPUT> duc_i;
        filterStageClass<l_INPUT> duc_q;

    public:
    
        //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
        void process(const std::array<DATA_T, l_INPUT>&  din_i,
                     std::array<DATA_T, l_OUTPUT>&  dout_i,
                     const std::array<DATA_T, l_INPUT>&  din_q,
                     std::array<DATA_T, l_OUTPUT>&  dout_q) {
        
        #pragma HLS INLINE
            duc_i.process(din_i, dout_i);
            duc_q.process(din_q, dout_q);
        }
}; 

#endif
