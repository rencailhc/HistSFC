//include SFC and PyramidT
#include "PyramidBuild.h"

void IOTbuild()
{
	CoordTrans pyramidtrans({ 19000,369000,-0.83,0 }, { 1,1,1,1 });
	double dimmax[] = { 20000, 370000, 25.68, 900 };
	double medians[] = { 19491.54, 369555.01,1.2,849.98 };
	//Pvalue_shift("E:/91_lod.csv", "E:/pyramid_91_3d.csv", pyramidtrans, dimmax, medians);	
	Pvalue_uni("E:/91_lod.csv", "E:/pyramid_91_3d_uni.csv", pyramidtrans, dimmax);
}