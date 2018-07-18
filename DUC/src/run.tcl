open_project -reset proj

add_files duc.cpp -cflags=-std=c++0x
add_files -tb duc_test.cpp -cflags=std=c++0x

set_top duc
#set_top filters
#set_top process_single_filter
open_solution -reset sol1

set_part {xc7vx690tffg1761-2}
create_clock -period 2 
set_clock_uncertainty 0

#set_part  {xc7k160tfbg484-1}
#create_clock -period 3
#set_part {xc7z020clg484-1}
#create_clock -period 5

#---- crucial directive for
# 1. to use fifo instead of BRAM
# 2. synthesis runs forever if not used
# depth can affect the II
#
# config_dataflow -default_channel fifo -fifo_depth 2
config_dataflow -default_channel fifo -fifo_depth 16

#csim_design
csynth_design
#cosim_design -rtl verilog -trace_level port
#cosim_design -rtl verilog -tool modelsim -trace_level port
export_design -evaluate verilog -format ip_catalog
#export_design -format sysgen
#export_design -evaluate verilog -format syn_dcp

exit

