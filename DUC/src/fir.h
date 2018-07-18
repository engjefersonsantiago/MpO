//______________________________________________________________________________
// fir.h:
// - various FIR module
//
//
// a. paek, sep 2013
//______________________________________________________________________________

#ifndef _FIR_H
#define _FIR_H

#include <cstdint>
#include <array>
#include <algorithm>
#include <iterator>

//______________________________________________________________________________
// filter precision
// - precisions used by all the FIR class
//______________________________________________________________________________

#define FIXED

#ifdef FIXED
#include <ap_fixed.h>

static constexpr auto NFRAC = 14;
typedef ap_fixed<3+NFRAC,  3, AP_TRN, AP_WRAP> DATA_T;
typedef ap_fixed<4+NFRAC,  4, AP_TRN, AP_WRAP> DATA2_T;
typedef ap_fixed<2+NFRAC,  2, AP_RND_CONV, AP_WRAP> COEF_T;
typedef ap_fixed<5+NFRAC,  5, AP_RND_CONV, AP_WRAP> PROD_T; // rounding helps ~20db dc offset
typedef ap_fixed<10+NFRAC, 10, AP_TRN, AP_WRAP> ACC_T;
#else
typedef float DATA_T;
typedef float DATA2_T;
typedef float COEF_T;
typedef float PROD_T;
typedef float ACC_T;
#endif

template<int l_INPUT, int l_OUTPUT, int l_TDL, int l_COEFF, int II_IntFactor>
class FirCommon {
    public:
        typedef std::array<DATA_T, l_INPUT> Din_arr_t;
        typedef std::array<DATA_T, l_OUTPUT> Dout_arr_t;
        typedef std::array<DATA_T, l_TDL> Sr_arr_t;
        typedef std::array<ACC_T, II_IntFactor> Acc_arr_t;
        typedef std::array<std::array<COEF_T, l_COEFF>, II_IntFactor> Coeff_arr_t;

        Sr_arr_t sr;
        Acc_arr_t acc;
        const Coeff_arr_t coeff;

		// SR: 1 input
        void shift_register(const DATA_T& din_, Sr_arr_t& sr_) {
        	//#pragma HLS
			std::copy(std::begin(sr_), std::end(sr_) - 1, std::begin(sr_) + 1);
			sr_[0] = din_;
        }

		// SR: Overload for 2 inputs
        void shift_register(const std::array<DATA_T, 2>& din_, Sr_arr_t& sr_) {
			std::copy(std::begin(sr_), std::end(sr_) - 1, std::begin(sr_) + 2);
			std::copy(std::begin(din_), std::end(din_), std::begin(sr_));
        }

        //constructor
        FirCommon(const Coeff_arr_t& cin) : coeff(cin) {
            #pragma HLS array_partition variable=sr  complete
			#pragma HLS array_partition variable=coeff complete
			#pragma HLS array_partition variable=acc complete
        }


};

struct MAC {
    // MAC engine: 1 operator
    ACC_T operator()(const DATA_T& din, const COEF_T& coef, const ACC_T&  acc ) {
    	PROD_T  prod   = din*coef;
        ACC_T   sum    = prod + acc;
        return sum;
    };

    // MAC engine overload: 2 operators for pre adder
    ACC_T operator()(const DATA_T&  din0, const DATA_T& din1, const COEF_T& coef, const ACC_T& acc) {
    	DATA2_T preadd = din0 + din1;
        PROD_T  prod   = preadd*coef;
        ACC_T   sum    = prod + acc;
        return sum;
    }
};

//______________________________________________________________________________
// Base FIR class
//
// - parameter:
//    -l_WHOLE:  number of taps/coeff
//    -l_SAMPLE: number of data samples to process
//    -II_GOAL:  initiation interval goal
//	  -D_FIR: derived class for static polymorphism
//______________________________________________________________________________

template<size_t l_WHOLE, size_t l_COEF, size_t l_SAMPLE, size_t II_GOAL, typename D_FIR>
class fir_base {

	protected:
		// Types
		typedef std::array<DATA_T, l_WHOLE> Data_arr_t;		// SR and data
		typedef std::array<DATA_T, l_SAMPLE> Sample_arr_t;	// sample data
		typedef std::array<COEF_T, l_COEF> Coeff_arr_t;	    // coeff

		// Constants
		static constexpr auto ODD    = l_WHOLE % 2;
		static constexpr auto l_HALF = l_WHOLE/2;

		// Data members
		Data_arr_t sr;
		ACC_T  acc;
		const Coeff_arr_t coeff;
        MAC Mac;

        // Init array function
        template<typename T, size_t N>
        static std::array<T, N> init_array(const T carr[N]) {
	        std::array<T, N> arr {};
		    std::copy(carr, &carr[N-1], arr.begin());
	        return arr;
        }

	public:

		//_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
		// filter
		void process (const DATA_T& din, DATA_T& dout) {
			#pragma HLS INLINE
			// using 'factor' instead of 'complete' uses BRAM instead of FF
			#pragma HLS array_reshape variable=sr  complete
			#pragma HLS array_reshape variable=coeff complete dim=0
			acc = 0;

			// Calls the specific method
			static_cast<D_FIR*>(this)->accum_impl(sr, acc, coeff);

			// SR
			std::copy(sr.begin(), sr.end() - 1, sr.begin() + 1);

			sr[0] = din;
			dout = acc;

		}

		//_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
		// filter frame
		//
		void process_frame(const Sample_arr_t& din, Sample_arr_t& dout)
		{
			// not much difference inline
			//#pragma HLS INLINE
			for (auto i = 0; i < din.size(); ++i) {
				#pragma HLS pipeline II=II_GOAL rewind
				process (din[i], dout[i]);
			}
		}

		//constructor
		fir_base(const Coeff_arr_t& cin) : coeff(cin) {}
		fir_base(const COEF_T cin[l_COEF]) : coeff(init_array<COEF_T, l_COEF>(cin)) {}

}; // fir_base

//______________________________________________________________________________
// single rate, non symmetric fir class
//
// - parameter:
//    -l_WHOLE:  number of taps/coeff
//    -l_SAMPLE: number of data samples to process
//    -II_GOAL:  initiation interval goal
//______________________________________________________________________________
template<size_t l_WHOLE, size_t l_COEF, size_t l_SAMPLE, size_t II_GOAL>
class nosym_fir : public fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, nosym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> {

	private:
		typedef fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, nosym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> Base;

	public:

		//_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
		// filter
		void accum_impl (const typename Base::Data_arr_t& sr, ACC_T& acc, const typename Base::Coeff_arr_t& coeff) {
			for (auto i = 0; i < coeff.size(); ++i) {
				acc = this->Mac(sr[i], coeff[i], acc);
			}
		}

		//constructor
		nosym_fir(const typename Base::Coeff_arr_t& cin) : Base(cin) {}
		nosym_fir(const COEF_T cin[l_COEF]) : Base(cin) {}

}; // nosym_fir

//______________________________________________________________________________
// single rate, symmetric fir classs
//
// - parameter:
//    -l_WHOLE:  number of taps
//    -l_SAMPLE: number of data sample
//    -II_GOAL:  initiation interval goal
//______________________________________________________________________________
template<size_t l_WHOLE, size_t l_COEF, size_t l_SAMPLE, size_t II_GOAL>
class sym_fir : public fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, sym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> {

	private:
		typedef fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, sym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> Base;

	public:

		void accum_impl (const typename Base::Data_arr_t& sr, ACC_T& acc, const typename Base::Coeff_arr_t& coeff) {

			for (auto i = 0; i < this->l_HALF; ++i) {
				acc = this->Mac(sr[i], sr[l_WHOLE-1-i], coeff[i], acc);
            }
			if (this->ODD) {
				acc = this->Mac(sr[this->l_HALF], coeff[this->l_HALF], acc);
            }
		}

		//constructor
		sym_fir(const typename Base::Coeff_arr_t& cin) : Base(cin) {}
		sym_fir(const COEF_T cin[l_COEF]) : Base(cin) {}

}; // sym_class

//______________________________________________________________________________
// single rate, halfband fir classs
//
// - parameter:
//    -l_WHOLE:  number of taps
//    -l_SAMPLE: number of data sample
//    -II_GOAL:  initiation interval goal
//______________________________________________________________________________
template<size_t l_WHOLE, size_t l_COEF, size_t l_SAMPLE, size_t II_GOAL>
class hb_fir : public fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, sym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> {

	private:
		typedef fir_base<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL, sym_fir<l_WHOLE, l_COEF, l_SAMPLE, II_GOAL>> Base;

	public:

		void accum_impl (const typename Base::Data_arr_t& sr, ACC_T& acc, const typename Base::Coeff_arr_t& coeff) {
			unsigned short k=0;

			for (int i = 0; i < l_COEF - 1; ++i) {
				acc = this->Mac(sr[k], sr[l_WHOLE-1-k], coeff[i], acc);
				k+=2;
			}

			// center tap
			acc = this->Mac(sr[this->l_HALF], 0, coeff[l_COEF-1], acc);
		}

		//constructor
		hb_fir(const typename Base::Coeff_arr_t& cin) : Base(cin) {}
		hb_fir(const COEF_T cin[l_COEF]) : Base(cin) {}

}; // hb_fir

template<int l_INPUT, int II_GOAL, int l_TDL, int l_COEFF, int II_IntFactor, typename D_FIR>
class FirInterp2Base {
	protected:
    	typedef FirCommon<l_INPUT, 2*l_INPUT, l_TDL, l_COEFF, II_IntFactor> Common;
    	Common firCommon_;

    public:

        void process(const DATA_T& din, std::array<DATA_T, 2>& dout) {
            #pragma HLS INLINE

            // Specific kernel
            auto derFir = static_cast<D_FIR*>(this);
            derFir->accum_impl(firCommon_.sr, firCommon_.acc, firCommon_.coeff);

            // SR
            firCommon_.shift_register(din, firCommon_.sr);

            // Ouput
            derFir->ouput_impl(dout, firCommon_.sr, firCommon_.acc, firCommon_.coeff);

        }

        void process_frame(const typename Common::Din_arr_t& din, typename Common::Dout_arr_t& dout) {
            std::array<DATA_T, 2> dout_tmp = {{ 0, 0 }};
            #pragma HLS array_partition variable=dout_tmp dim=1
            
            for (auto i = 0; i < din.size(); ++i) {
                #pragma HLS pipeline II=II_GOAL rewind
                process(din[i], dout_tmp);
                std::copy(dout_tmp.begin(), dout_tmp.end(), &dout[2*i]);
            }
        }

        //constructor
        FirInterp2Base(const typename Common::Coeff_arr_t& cin) : firCommon_(cin) {}
};

//______________________________________________________________________________
// interpolate by 2 FIR class
//
// - make sure the coeff is ODD length, so both subfilter will be symmetric
// - parameter:
//    -l_WHOLE:   number of taps of prototype filter (not the decomposed subfilter)
//    -l_INPUT:  number of data samples
//    -II_GOAL:   initiation interval goal
//
//______________________________________________________________________________
template<int l_WHOLE> struct FirInterp2Cfg {
    static constexpr auto INTERP_FACTOR = 2;
    static constexpr auto L_SUB         = l_WHOLE/INTERP_FACTOR;  // 32
    static constexpr auto ODD           = l_WHOLE % 2;            // 1
    static constexpr auto l_TDL         = L_SUB + ODD;            // 33
    static constexpr auto l_NONZERO     = L_SUB/2 + ODD;          // 17
};

template<int l_WHOLE, int l_INPUT, int II_GOAL>
class FirInterp2 : public FirInterp2Base<l_INPUT, II_GOAL,
		FirInterp2Cfg<l_WHOLE>::l_TDL,
		FirInterp2Cfg<l_WHOLE>::l_NONZERO,
		FirInterp2Cfg<l_WHOLE>::INTERP_FACTOR,
		FirInterp2<l_WHOLE, l_INPUT, II_GOAL>> {

    protected:

        typedef FirInterp2Cfg<l_WHOLE> Config;
        typedef FirInterp2Base<l_INPUT, II_GOAL, Config::l_TDL,
            Config::l_NONZERO, Config::INTERP_FACTOR,
			FirInterp2<l_WHOLE, l_INPUT, II_GOAL>> Base;
        typedef typename Base::Common Common;

    public:

        void accum_impl(const typename Common::Sr_arr_t& sr,
                        typename Common::Acc_arr_t& acc,
                        const typename Common::Coeff_arr_t& coeff) {
            acc = {{ 0, 0 }};
            MAC mac;
            for (auto i = 0; i < Config::l_NONZERO - Config::ODD; ++i) {
                // even number of taps, has one more than ODD one
                acc[0] = mac(sr[i], sr[Config::l_TDL-1-i], coeff[0][i], acc[0]);

                // ODD number of taps
                acc[1] = mac(sr[i], sr[Config::l_TDL-1-Config::ODD-i], coeff[1][i], acc[1]);
            }

            // center tap
            acc[0] = mac(sr[Config::l_NONZERO-1], coeff[0][Config::l_NONZERO-1], acc[0]);
        }

        void ouput_impl(std::array<DATA_T, 2>& dout,
        				const typename Common::Sr_arr_t&,
                        const typename Common::Acc_arr_t& acc,
                        const typename Common::Coeff_arr_t&)  {
        	dout = {{ acc[0], acc[1] }};
        }

        // initialize coefficient for polyphase decomposition
        // - test out even length
        const typename Common::Coeff_arr_t init (const COEF_T coef_in[l_WHOLE]) const {
            static constexpr auto gain = 2;  // only for DUC
            typename Common::Coeff_arr_t coeff;
            for (int i = 0; i < Config::l_NONZERO - Config::ODD; ++i) {
                for (int k = 0; k < Config::INTERP_FACTOR; k++) {
                    coeff[k][i] = gain*coef_in[2*i+k];
                }
            }
            // number of taps is one greater for ODD phase filter
            if (Config::ODD) {
                static constexpr auto i = Config::l_NONZERO - 1;
                coeff[0][i] = gain * coef_in[2*i];
            }
            return coeff;
        }

        //constructor
        FirInterp2(const COEF_T cin[l_WHOLE]) : Base(init(cin)) {}

}; // interp2_class

//______________________________________________________________________________
// interpolate by 2 halfband fir class
//
// - parameter:
//    -l_WHOLE:   number of taps
//    -l_INPUT:  number of input data samples
//    -II_GOAL:   initiation interval goal
//
//
//    -l_TDL:     number of stages in tap delay line
//    -l_NONZERO: number of unique coefficients
//______________________________________________________________________________
template<int l_WHOLE> struct FirInterp2HbCfg {
    static constexpr auto l_TDL     = (l_WHOLE+1)/2; // (23+1)/2 = 12
    static constexpr auto l_NONZERO = 1+l_TDL/2;     // 1 + 12/2 = 7
};

template<int l_WHOLE, int l_INPUT, int II_GOAL>
class FirInterp2Hb : public FirInterp2Base<l_INPUT, II_GOAL,
		FirInterp2HbCfg<l_WHOLE>::l_TDL,
		FirInterp2HbCfg<l_WHOLE>::l_NONZERO, 1,
		FirInterp2Hb<l_WHOLE, l_INPUT, II_GOAL>> {

    protected:
		typedef FirInterp2HbCfg<l_WHOLE> Config;
		typedef FirInterp2Base<l_INPUT, II_GOAL, Config::l_TDL,
			Config::l_NONZERO, 1,
			FirInterp2Hb<l_WHOLE, l_INPUT, II_GOAL>> Base;
		typedef typename Base::Common Common;

    public:

        //_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        // filter
        void accum_impl(const typename Common::Sr_arr_t& sr,
                        typename Common::Acc_arr_t& acc,
                        const typename Common::Coeff_arr_t& coeff) {
            acc = {{ 0 }};
            MAC mac;
            for (auto i = 0; i < Config::l_NONZERO - 1; ++i) {
                acc[0] = mac(sr[i], sr[Config::l_TDL-1-i], coeff[0][i], acc[0]);
            }

        }

        void ouput_impl(std::array<DATA_T, 2>& dout,
        				const typename Common::Sr_arr_t& sr,
                        const typename Common::Acc_arr_t& acc,
                        const typename Common::Coeff_arr_t& coeff) {
            dout = {{ acc[0], coeff[0][Config::l_NONZERO-1]*sr[Config::l_NONZERO-1] }};
        }

        //_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        // initialize coefficient for halfband filter
        const typename Common::Coeff_arr_t init(const COEF_T cin[]) const {
            static constexpr auto gain = 2;  // only for DUC
            typename Common::Coeff_arr_t coeff;
            for (auto i = 0; i < Config::l_NONZERO; ++i) {
                coeff[0][i]  = gain*((i == Config::l_NONZERO - 1) ? cin[2*i-1] : cin[2*i]);
            }
            return coeff;
        }

        //constructor
        FirInterp2Hb(const COEF_T cin[l_WHOLE]) : Base(init(cin)) {}

}; // interp2_hb_class

//______________________________________________________________________________
// decimate by 2, symmetric fir classs
//
// - parameter:
//    -l_WHOLE:   number of taps
//    -l_OUTPUT:  number of output data samples
//    -II_GOAL:   initiation interval goal
//______________________________________________________________________________

template<int l_WHOLE, int l_OUTPUT, int II_GOAL>
class decim2_class {

    private:
        static constexpr auto ODD    = l_WHOLE % 2;
        static constexpr auto l_HALF = l_WHOLE/2;
        static constexpr auto l_COEF = l_WHOLE/2 + ODD;

        DATA_T sr[l_WHOLE];
        ACC_T  acc;
        COEF_T coeff[l_COEF];

    public:

        //_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        // filter
        void process( DATA_T din0, DATA_T din1, DATA_T* dout) {

            // using 'factor' instead of 'complete' uses BRAM instead of FF
            #pragma HLS array_reshape variable=sr  complete
            #pragma HLS array_reshape variable=coeff complete dim=0

            acc = 0;

            for (int i=0; i<l_HALF; i++)
                acc = MAC(sr[i], sr[l_WHOLE-1-i], coeff[i], acc);

            if (ODD)
                acc = MAC(sr[l_HALF], coeff[l_HALF], acc);

            for (int i=l_WHOLE-1; i>1 ; i--)
                sr[i] = sr[i-2];

            sr[1] = din0;
            sr[0] = din1; // most recent sample

            *dout = acc;

        }

        //_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
        // filter frame
        void process_frame(DATA_T din[2*l_OUTPUT], DATA_T dout[l_OUTPUT]) {

            DATA_T din0, din1;

            for (int i=0; i<l_OUTPUT; i++ ) {
                #pragma HLS pipeline II=II_GOAL rewind
                din0 = din[2*i];
                din1 = din[2*i+1];
                process (din0, din1, dout[i]);
            }
        }

        //_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        // initialize coefficient
        void init(const COEF_T cin[l_WHOLE] ) {

            for (int i=0; i<l_COEF-1; i++) {
                coeff[i]  = cin[i];
            }
        }

        //constructor
        decim2_class(const COEF_T cin[l_WHOLE]) {
            init(cin);
        }

        decim2_class(void) {}

}; // decim2_class

//______________________________________________________________________________
// decimate by 2, halfband fir classs
//
// - parameter:
//    -l_WHOLE:   number of taps
//
//______________________________________________________________________________

template<int l_WHOLE>
class hb_decim2_class {

    protected:
        static const int l_HALF = l_WHOLE/2;
        static const int l_COEF = (l_WHOLE+1)/4+1; // number of nonzero coeff

        DATA_T sr[l_WHOLE];
        ACC_T  acc;

    public:

        //_  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _  _
        // filter

        void process ( DATA_T din0, DATA_T din1, DATA_T* dout, const COEF_T cin[l_COEF]) {

            // using 'factor' instead of 'complete' uses BRAM instead of FF
            #pragma HLS array_reshape variable=sr  complete
            #pragma HLS array_reshape variable=cin complete dim=0

            acc = 0;
            unsigned short k = 0;

            for (int i=0; i<l_COEF-1; i++) {
                acc = this->MAC_preadd (sr[k], sr[l_WHOLE-1-k], cin[i], acc);
                k+=2;  // skipping zero coefficient
            }

            // center tap
            acc = this->MAC_preadd (sr[l_HALF], 0, cin[l_COEF-1], acc);

            for (int i=l_WHOLE-1; i>1 ; i--)
                sr[i] = sr[i-2];

            sr[1] = din0;
            sr[0] = din1; // most recent sample

            *dout = acc;

        }

}; // hb_decim2_class

#endif
