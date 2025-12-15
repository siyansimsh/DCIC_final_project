# Temporary placeholder constraints for simulation/analysis only
# Update with real board pinout (LOC) and actual I/O standards before tapeout/board bring-up.

# Clock constraint (adjust period and port name if different)
create_clock -name sys_clk -period 10.000 [get_ports clk]

# Generic I/O standard to suppress NSTD (change to real voltage/standard later)
set_property IOSTANDARD LVCMOS18 [get_ports *]

# Placeholder input delays (relative to sys_clk)
# Replace 2.0/0.5 with real tCO + board skew numbers per interface.
set input_ports [list \
    rst_n \
    in_valid \
    adc1_re adc1_im \
    adc2_re adc2_im \
    h00_re h00_im \
    h01_re h01_im \
    h10_re h10_im \
    h11_re h11_im]
set_input_delay -clock [get_clocks sys_clk] -max 2.0 [get_ports $input_ports]
set_input_delay -clock [get_clocks sys_clk] -min 0.5 [get_ports $input_ports]

# Placeholder output delays (relative to sys_clk)
# Replace 2.0/-0.5 with sink tsu/th + board skew.
set output_ports [list \
    out_valid \
    rx_bits[*] \
    zf_valid_tap \
    cpe_phase_valid_tap \
    zf_x1_re_tap zf_x1_im_tap \
    zf_x2_re_tap zf_x2_im_tap \
    cpe_phase_err_tap \
    nco_acc_tap]
set_output_delay -clock [get_clocks sys_clk] -max 2.0 [get_ports $output_ports]
set_output_delay -clock [get_clocks sys_clk] -min -0.5 [get_ports $output_ports]

# Optional: ignore timing on async reset if treated as asynchronous
#set_false_path -from [get_ports rst_n]
