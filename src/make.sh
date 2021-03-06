#!/bin/bash

ulimit -s unlimited
#ulimit -s 102400

cd ..
make
cd ./src

read

#./hydrocode.out Sod_10_test2 Sod_10_test2/Sod_10_test_ROE 2 1_Roe t1

#./hydrocode.out Sod_test Sod_test/Sod_test_ROE	1 1_Roe		  free
#./hydrocode.out Sod_test Sod_test/Sod_test	1 1_Riemann_exact free
#./hydrocode.out Sod_test Sod_test/Sod_test	1 2_GRP		  free
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test_ROE 2 1_Roe		  Sod
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test	 2 1_Riemann_exact Sod
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test	 2 2_GRP		  Sod
#./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test	 2 2_GRP_2D		  Sod

#./hydrocode.out odd_even odd_even/odd_even_Roe 2 1_Roe 	  odd_even
#./hydrocode.out odd_even odd_even/odd_even	2 1_Riemann_exact odd_even

#./hydrocode.out odd_even odd_even/odd_even	2 1_Riemann_exact Sod

#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_O#./hydrocode.out 2DSI_1281	2DSI_1281/2DSI_1281	2 1_Riemann_exact RMIBLIQUE_Roe	2 1_Roe			oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 1_Riemann_exact	oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE	2 2_GRP			oblique_periodic
#./hydrocode.out NEW_TEST_OBLIQUE_SMALL NEW_TEST_OBLIQUE_SMALL/NEW_TEST_OBLIQUE_SMALL	2 1_Riemann_exact			oblique_periodic

#./hydrocode.out NEW_TEST_Sod NEW_TEST_Sod/NEW_TEST_Sod	2 1_Riemann_exact	oblique_periodic
#./hydrocode.out NEW_TEST_Shock NEW_TEST_Shock/NEW_TEST_Shock	2 1_Riemann_exact	oblique_periodic

#./hydrocode.out NEW_TEST NEW_TEST/NEW_TEST_Roe	2 1_Roe			free
#./hydrocode.out NEW_TEST NEW_TEST/NEW_TEST	2 1_Riemann_exact	free
#./hydrocode.out NEW_TEST_BIG NEW_TEST_BIG/NEW_TEST_BIG_Roe	2 1_Roe			free

#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_Roe	2 1_Roe			free
#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad		2 1_Riemann_exact	free
#./hydrocode.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad		2 2_GRP			free
#./hydrocode.out RP2D_shear_Quad	RP2D_shear_Quad/RP2D_shear_Quad		2 2_GRP			RMI_S
#./hydrocode.out RP2D_shear_Quad	RP2D_shear_Quad/RP2D_shear_Quad_G2D		2 2_GRP_2D		RMI_S
#./hydrocode.out RP2D_shear_Quad_6	RP2D_shear_Quad_6/RP2D_shear_Quad_GRP	2 2_GRP		RMI_S
#./hydrocode.out RP2D_shear_Quad_6	RP2D_shear_Quad_6/RP2D_shear_Quad_G2D	2 2_GRP_2D		RMI_S
#./hydrocode.out RP2D_shear_Quad_6700	RP2D_shear_Quad_6700/RP2D_shear_Quad_G2D	2 2_GRP_2D		RMI_S

#./hydrocode.out RMI/RMI_81	RMI/RMI_81	2 1_Riemann_exact RMI
#./hydrocode.out RMI/RMI_81	RMI/RMI_81	2 2_GRP		  RMI
#./hydrocode.out RMI/RMI_321	RMI/RMI_321	2 1_Riemann_exact RMI
#./hydrocode.out RMI/RMI_321	RMI/RMI_321	2 2_GRP		  RMI
#./hydrocode.out RMI_256	RMI_256	2 1_Riemann_exact RMI
#./hydrocode.out RMI_256	RMI_256_3	2 2_GRP		  RMI
#./hydrocode.out RMI_128	RMI_128_3	2 2_GRP		  RMI
#./hydrocode.out RMI_192	RMI_192_3	2 2_GRP		  RMI_S
#./hydrocode.out RMI_384	RMI_384_3	2 2_GRP		  RMI_S
#./hydrocode.out RMI_512	RMI_512_1.5	2 2_GRP		  RMI
#./hydrocode.out RMI_768	RMI_768_1.5	2 2_GRP		  RMI_S
#./hydrocode.out RMI_easy/128	RMI_easy/128_GRP	2 1_Riemann_exact RMI
#./hydrocode.out RMI_easy/128	RMI_easy/128_GRP	2 2_GRP		  RMI
#./hydrocode.out RMI_easy/128_small	RMI_easy/128_small	2 2_GRP		  RMI
#./hydrocode.out RMI_easy/128	RMI_easy/128_G2D	2 2_GRP_2D		  RMI
#./hydrocode.out RMI_easy/128_small	RMI_easy/128_small_G2D	2 2_GRP_2D		  RMI
#./hydrocode.out RMI_easy/128_RHOdiff	RMI_easy/128_RHOdiff_GRP	2 2_GRP		  RMI
#./hydrocode.out RMI_easy/128_RHOdiff	RMI_easy/128_RHOdiff_G2D	2 2_GRP_2D		  RMI

#./hydrocode.out 2DSI/2DSI_1281	2DSI_1281/2DSI_1281	2 1_Riemann_exact RMI
#nohup ./hydrocode.out 2DSI/2DSI_A0K1	2DSI/2DSI_A0K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A8K1	2DSI/2DSI_A8K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A16K1	2DSI/2DSI_A16K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A20K1	2DSI/2DSI_A20K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A40K1	2DSI/2DSI_A40K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A80K1	2DSI/2DSI_A80K1_fix	2 1_Riemann_exact RMI &
#nohup ./hydrocode.out 2DSI/2DSI_A200K1 2DSI/2DSI_A200K1_fix	2 1_Riemann_exact RMI &

#./hydrocode.out 2D_Shock–interface/2DSI_A16K1	2D_Shock–interface/2DSI_A16K1	2 2_GRP RMI

#./hydrocode.out A3_shell	A3_shell/God1	2 1_Riemann_exact Shell
#./hydrocode.out A3_shell	A3_shell/GRP2	2 2_GRP Shell
#./hydrocode.out A3_shell	A3_shell/GRP2_G2D	2 2_GRP_2D Shell
#./hydrocode.out A3_shell_4	A3_shell_4/God1	2 1_Riemann_exact free
#./hydrocode.out A3_shell_4	A3_shell_4/GRP2	2 2_GRP RMI_S
#./hydrocode.out A3_shell_4	A3_shell_4/GRP2_G2D	2 2_GRP_2D RMI_S
#./hydrocode.out Bubble	Bubble/God1	2 1_Riemann_exact Sod
#./hydrocode.out Bubble	Bubble/GRP2	2 2_GRP Sod
#./hydrocode.out Bubble_He	Bubble_He/God1	2 1_Riemann_exact Sod
#./hydrocode.out Bubble_He	Bubble_He/GRP2	2 2_GRP Sod
#./hydrocode.out Bubble_He_RHOdiff	Bubble_He_RHOdiff/GRP2	2 2_GRP Sod
#./hydrocode.out Bubble_R22	Bubble_R22/God1	2 1_Riemann_exact Sod
#./hydrocode.out Bubble_R22	Bubble_R22/GRP2	2 2_GRP Sod
#./hydrocode.out Bubble_R22_RHOdiff	Bubble_R22_RHOdiff/GRP2		2 2_GRP	Sod
#./hydrocode.out Bubble_R22_RHOdiff	Bubble_R22_RHOdiff/GRP_G2D	2 2_GRP_2D Sod
#./hydrocode.out Bubble_He_6	Bubble_He_6/GRP2	2 2_GRP Sod
#./hydrocode.out Bubble_SF6/b	Bubble_SF6/b_185/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/c	Bubble_SF6/c_185/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/d	Bubble_SF6/d_185/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/e	Bubble_SF6/e_185/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/OLD_800/b	Bubble_SF6/b_800_18N/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/OLD_800/c	Bubble_SF6/c_800_18N/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/OLD_800/d	Bubble_SF6/d_800_18N/GRP2	2 2_GRP Shock_Bubble
#./hydrocode.out Bubble_SF6/OLD_800/e	Bubble_SF6/e_800_18N/GRP2	2 2_GRP Shock_Bubble

#./hydrocode.out Bubble_He		Bubble_He/GRP	2 2_GRP Sod
#./hydrocode.out Bubble_R22	Bubble_R22/GRP	2 2_GRP Sod
#./hydrocode.out Bubble_He		Bubble_He/GRP2D	2 2_GRP_2D Sod
#./hydrocode.out Bubble_R22	Bubble_R22/GRP2D	2 2_GRP_2D Sod

#./hydrocode.out Bubble_He_X4		Bubble_He_X4/GRP	2 2_GRP Sod
#./hydrocode.out Bubble_R22_X4	Bubble_R22_X4/GRP	2 2_GRP Sod
#./hydrocode.out Bubble_He_X4		Bubble_He_X4/GRP2D	2 2_GRP_2D Sod
#./hydrocode.out Bubble_R22_X4	Bubble_R22_X4/GRP2D	2 2_GRP_2D Sod

#./hydrocode.out Rare_isolate	Rare_isolate/God1	2 1_Riemann_exact free
#./hydrocode.out Rare_isolate	Rare_isolate/God1	2 2_GRP oblique_periodic
#./hydrocode.out Tangent_Discon_OBLIQUE	Tangent_Discon_OBLIQUE/Tangent_Discon_OBLIQUE_nofix	2 1_Riemann_exact Shear
#./hydrocode.out 1D_Shock-interface_5	1D_Shock-interface_5/God1	2 1_Riemann_exact free
#./hydrocode.out 1D_Shock-interface_5	1D_Shock-interface_5/GRP2	2 2_GRP oblique_periodic
#./hydrocode.out 1D_Shock-interface_5_NEW	1D_Shock-interface_5_NEW/God1	2 1_Riemann_exact free
#./hydrocode.out 1D_Shock-interface_5_NEW	1D_Shock-interface_5_NEW/GRP2	2 2_GRP oblique_periodic
#./hydrocode.out 1D_Shock-interface_5_NEW3	1D_Shock-interface_5_NEW3/God1	2 1_Riemann_exact free
#./hydrocode.out 2M_interface	2M_interface/2M_interface	2 1_Riemann_exact free
#./hydrocode.out 2M_interface	2M_interface/2M_interface	2 2_GRP oblique_periodic
#./hydrocode.out 2M_interface_ST	2M_interface_ST/2M_interface_ST	2 2_GRP free
#./hydrocode.out 2M_diff_P	2M_diff_P/2M_diff_P	2 1_Riemann_exact free
#./hydrocode.out Numerical_Failures	Numerical_Failures/Numerical_Failures	2 2_GRP free

#./hydrocode.out Is_Vortex/u_0		Is_Vortex/u_0/GRP2		2 2_GRP	Vortex
#./hydrocode.out Is_Vortex/u_0_5	Is_Vortex/u_0_5/GRP2		2 2_GRP	Vortex
#./hydrocode.out Is_Vortex/u_0_5_v_0_5	Is_Vortex/u_0_5_v_0_5/GRP2	2 2_GRP	Vortex
#./hydrocode.out Is_Vortex/u_0		Is_Vortex/u_0/GRP2_G2D	2 2_GRP_2D	Vortex
#./hydrocode.out Is_Vortex/u_0_5	Is_Vortex/u_0_5/GRP2_G2D 	2 2_GRP_2D	Vortex
#./hydrocode.out Is_Vortex/u_0_5_v_0_5	Is_Vortex/u_0_5_v_0_5/GRP2_G2D 	2 2_GRP_2D	Vortex

#./hydrocode.out Accuracy_test/20		Accuracy_test/20/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/40		Accuracy_test/40/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/80		Accuracy_test/80/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/160		Accuracy_test/160/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/320		Accuracy_test/320/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/640		Accuracy_test/640/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/1280		Accuracy_test/1280/GRP2		2 2_GRP	Vortex
#./hydrocode.out Accuracy_test/40		Accuracy_test/40/GRPQ1D		2 2_GRP_2D	Vortex
#./hydrocode.out Accuracy_test/80		Accuracy_test/80/GRPQ1D		2 2_GRP_2D	Vortex
#./hydrocode.out Accuracy_test/160		Accuracy_test/160/GRPQ1D	2 2_GRP_2D	Vortex
#./hydrocode.out Accuracy_test/320		Accuracy_test/320/GRPQ1D	2 2_GRP_2D	Vortex
#./hydrocode.out Accuracy_test/640		Accuracy_test/640/GRPQ1D	2 2_GRP_2D	Vortex


#./hydrocode.out RMI_easy/128_small_new		RMI_easy/128_small_new_G2D	2 2_GRP_2D		  RMI
#./hydrocode.out RMI_easy/128_small_new		RMI_easy/128_small_new_P1D	2 2_GRP		  RMI

gprof -b -A -p -q hydrocode.out  gmon.out > pg
