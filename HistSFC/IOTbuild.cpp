/*Building the Pyramid-Technique using Oracle IOT*/

#include "PyramidBuild.h"

void IOTbuild()
{
	CoordTrans pyramidtrans({ 19000,369000,-0.83,0 }, { 1,1,1,1 });
	double dimmax[] = { 20000, 370000, 25.68, 900 };
	double medians[] = { 19491.54, 369555.01,1.2,849.98 };
	//Pvalue_shift("E:/91_lod.csv", "E:/pyramid_91_3d.csv", pyramidtrans, dimmax, medians);	
	//Pvalue_uni("E:/91_lod.csv", "E:/pyramid_91_3d_uni.csv", pyramidtrans, dimmax);
	CoordTrans sim1trans({ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }, { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 });
	double sim1dimmax[] = { 4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096 };
	double sim1midians[] = { 1757,2186,2606,2253,2074,2106,1383,1753,1436, 2910, 2452, 1217, 2685,2046, 2180, 3089 };
	//Pvalue_uni("E:/sigmod/simres/pt/sim1_data.csv", "E:/sigmod/simres/pt/pt_sim1_uni_data.csv",sim1trans,sim1dimmax);
	//Pvalue_shift("E:/sigmod/simres/pt/sim1_data.csv", "E:/sigmod/simres/pt/pt_sim1_shift_data.csv", sim1trans, sim1dimmax,sim1midians);

	CoordTrans sim2trans({ 0,0,0,0,0,0 }, { 1,1,1,1,1,1 });
	double sim2dimmax[] = { 1 << 20,1 << 20,1 << 20,1 << 20,1 << 20,1 << 20 };
	double sim2midians_ind[] = { 524555,524328.5,176956.5,272217,165051.5,879803 };
	double sim2midians_cor[] = { 524301,527249.5,176467,269097,164903,880241 };
	Pvalue_uni("E:/sigmod/simres/sim2/realsim_cor.csv", "E:/sigmod/simres/pt/pt_sim2_uni_cor_data.csv", sim2trans, sim2dimmax);
	Pvalue_shift("E:/sigmod/simres/sim2/realsim_cor.csv", "E:/sigmod/simres/pt/pt_sim2_shift_cor_data.csv", sim2trans, sim2dimmax, sim2midians_cor);
}
