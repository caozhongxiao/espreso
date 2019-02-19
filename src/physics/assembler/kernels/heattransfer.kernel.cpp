
#include "heattransfer.kernel.h"

#include "basis/evaluator/evaluator.h"
#include "esinfo/eslog.hpp"
#include "config/ecf/physics/heattransfer.h"

#include <cmath>
#include <algorithm>

using namespace espreso;

static void convectionMaterialParameters(
		const ConvectionConfiguration &convection,
		esint csize, double *coordinates, double time, double temp, double extTemp,
		double &rho, double &dynamicViscosity, double &dynamicViscosityTemp, double &heatCapacity, double &thermalConductivity)
{
	double  gas_constant;

	switch (convection.fluid) {
	case ConvectionConfiguration::FLUID::AIR: {
		gas_constant = 286.9;
		convection.absolute_pressure.evaluator->evalVector(1, csize, coordinates, &extTemp, time, &rho);
		rho /= (gas_constant * extTemp);

		if ((extTemp >=200) && (extTemp <= 1600)){
			heatCapacity = 1047.63657-0.372589265*extTemp+9.45304214E-4*pow(extTemp,2.0)-6.02409443E-7*pow(extTemp,3.0)+1.2858961E-10*pow(extTemp,4.0);
		}else if (extTemp < 200){
			heatCapacity = 1047.63657-0.372589265*200.0+9.45304214E-4*pow(200.0,2.0)-6.02409443E-7*pow(200.0,3.0)+1.2858961E-10*pow(200.0,4.0);
		}else if (extTemp > 1600){
			heatCapacity = 1047.63657-0.372589265*1600.0+9.45304214E-4*pow(1600.0,2.0)-6.02409443E-7*pow(1600.0,3.0)+1.2858961E-10*pow(1600.0,4.0);
		}

		if ((extTemp >=200) && (extTemp <= 1600)){
			thermalConductivity = -0.00227583562+1.15480022E-4*extTemp-7.90252856E-8*pow(extTemp,2.0)+4.11702505E-11*pow(extTemp,3.0)-7.43864331E-15*pow(extTemp,4.0);
		}else if (extTemp < 200){
			thermalConductivity = -0.00227583562+1.15480022E-4*200.0-7.90252856E-8*pow(200.0,2.0)+4.11702505E-11*pow(200.0,3.0)-7.43864331E-15*pow(200.0,4.0);
		}else if (extTemp > 1600){
			thermalConductivity =  -0.00227583562+1.15480022E-4*1600.0-7.90252856E-8*pow(1600.0,2.0)+4.11702505E-11*pow(1600.0,3.0)-7.43864331E-15*pow(1600.0,4.0);
		}

		if ((extTemp >=200) && (extTemp <= 1600)){
		dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * extTemp - 7.69429583E-11 * pow(extTemp,2.0) + 4.6437266E-14 * pow(extTemp,3.0) - 1.06585607E-17 * pow(extTemp,4.0);
		}else if (extTemp < 200){
			dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (extTemp > 1600){
			dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}


		if ((temp >=200) && (temp <= 1600)){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * temp - 7.69429583E-11 * pow(temp,2.0) + 4.6437266E-14 * pow(temp,3.0) - 1.06585607E-17 * pow(temp,4.0);
		}else if (temp < 200){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (temp > 1600){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}

	}break;
	case ConvectionConfiguration::FLUID::WATER:{

		if ((extTemp >=273) && (extTemp <= 283)){
			rho = 972.7584 + 0.2084 *extTemp - 4.0E-4 * pow(extTemp,2.0);
		}else if((extTemp >283) && (extTemp <= 373)){
			rho = 345.28 + 5.749816 * extTemp - 0.0157244 * pow(extTemp,2.0) + 1.264375E-5 * pow(extTemp,3.0);
		}else if (extTemp < 273){
			rho = 972.7584 + 0.2084 *273.0 - 4.0E-4 * pow(273.0,2.0);
		}else if (extTemp > 373){
			rho = 345.28 + 5.749816 * 373.0 - 0.0157244 * pow(373.0,2.0) + 1.264375E-5 * pow(373.0,3.0);
		}


		if ((extTemp >= 265) && (extTemp <= 293)){
			dynamicViscosity = 5.948859 - 0.08236196 * extTemp + 4.287142E-4 * pow(extTemp,2.0) - 9.938045E-7 * pow(extTemp,3.0) + 8.65316E-10 * pow(extTemp,4.0);
		}else if((extTemp >293) && (extTemp <= 353)){
			dynamicViscosity = 	0.410191 - 0.004753985 * extTemp + 2.079795E-5 * pow(extTemp,2.0) - 4.061698E-8 *  pow(extTemp,3.0) + 2.983925E-11 * pow(extTemp,4.0);
		}else if((extTemp >353) && (extTemp <= 423)){
			dynamicViscosity = 0.03625638 - 3.265463E-4 * extTemp + 1.127139E-6 * pow(extTemp,2.0) - 1.75363E-9 * pow(extTemp,3.0) + 1.033976E-12 * pow(extTemp,4.0);
		}else if (extTemp < 265){
			dynamicViscosity = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (extTemp > 423){
			dynamicViscosity = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}


		if ((extTemp >= 275) && (extTemp <= 370)){
			thermalConductivity = -0.9003748 + 0.008387698 * extTemp - 1.118205E-5 * pow(extTemp,2.0);
		}else if (extTemp < 275){
			thermalConductivity = -0.9003748 + 0.008387698 * 275.0 - 1.118205E-5 * pow(275.0,2.0);
		}else if (extTemp > 370){
			thermalConductivity = -0.9003748 + 0.008387698 * 370.0 - 1.118205E-5 * pow(370.0,2.0);
		}

		if ((extTemp >= 293) && (extTemp <= 373)){
			heatCapacity = 4035.841 + 0.492312 * extTemp;
		}else if (extTemp < 293){
			heatCapacity = 4035.841 + 0.492312 * 293.0;
		}else if (extTemp > 373){
			heatCapacity = 4035.841 + 0.492312 * 373.0;
		}


		if ((temp >= 265) && (temp <= 293)){
			dynamicViscosityTemp = 5.948859 - 0.08236196 * temp + 4.287142E-4 * pow(temp,2.0) - 9.938045E-7 * pow(temp,3.0) + 8.65316E-10 * pow(temp,4.0);
		}else if((temp >293) && (temp <= 353)){
			dynamicViscosityTemp = 	0.410191 - 0.004753985 * temp + 2.079795E-5 * pow(temp,2.0) - 4.061698E-8 *  pow(temp,3.0) + 2.983925E-11 * pow(temp,4.0);
		}else if((temp >353) && (temp <= 423)){
			dynamicViscosityTemp = 0.03625638 - 3.265463E-4 * temp + 1.127139E-6 * pow(temp,2.0) - 1.75363E-9 * pow(temp,3.0) + 1.033976E-12 * pow(temp,4.0);
		}else if (temp < 265){
			dynamicViscosityTemp = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (temp > 423){
			dynamicViscosityTemp = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}

	}break;
	case ConvectionConfiguration::FLUID::ENGINE_OIL:{


		if ((extTemp >=273) && (extTemp <= 433)){
			rho = 1068.70404 - 0.6393421 * extTemp + 7.34307359E-5 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			rho = 1068.70404 - 0.6393421 * 273.0 + 7.34307359E-5 * pow(273.0,2.0);
		}else if (extTemp > 433){
			rho = 1068.70404 - 0.6393421 * 433.0 + 7.34307359E-5 * pow(433.0,2.0);
		}

		if ((extTemp >= 273) && (extTemp <= 433)){
			thermalConductivity = 0.192223542 - 2.0637987E-4 * extTemp + 1.54220779E-7 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			thermalConductivity = 0.192223542 - 2.0637987E-4 * 273.0 + 1.54220779E-7 * pow(273.0,2.0);
		}else if (extTemp > 433){
			thermalConductivity =0.192223542 - 2.0637987E-4 * 433.0 + 1.54220779E-7 * pow(433.0,2.0);
		}


		if ((extTemp >= 273) && (extTemp <= 433)){
			heatCapacity = 761.405625 + 3.47685606 * extTemp + 0.00115530303 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			heatCapacity = 761.405625 + 3.47685606 * 273.0 + 0.00115530303 * pow(273.0,2.0);
		}else if (extTemp > 433){
			heatCapacity = 761.405625 + 3.47685606 * 433.0 + 0.00115530303 * pow(433.0,2.0);
		}


		if ((extTemp >= 273) && (extTemp <= 353)){
			dynamicViscosity = 42669.28688622 - 741.1718801282 * extTemp + 5.360521287088 * pow(extTemp,2.0) - 0.02066027676164 * pow(extTemp,3.0) + 4.47491538052E-5 * pow(extTemp,4.0) - 5.164053479202E-8 * pow(extTemp,5.0) + 2.48033770504E-11 * pow(extTemp,6.0);
		}else if ((extTemp > 353) && (extTemp <= 433 )){
			dynamicViscosity = 4.94593941 - 0.0351869631 * extTemp + 8.37935977E-5 * pow(extTemp,2.0) - 6.67125E-8 * pow(extTemp,3.0);

		}else if (extTemp < 273){
			dynamicViscosity = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (extTemp > 433){
			dynamicViscosity = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}

		if ((temp >= 273) && (temp <= 353)){
			dynamicViscosityTemp = 42669.28688622 - 741.1718801282 * temp + 5.360521287088 * pow(temp,2.0) - 0.02066027676164 * pow(temp,3.0) + 4.47491538052E-5 * pow(temp,4.0) - 5.164053479202E-8 * pow(temp,5.0) + 2.48033770504E-11 * pow(temp,6.0);
		}else if ((temp > 353) && (temp <= 433 )){
			dynamicViscosityTemp = 4.94593941 - 0.0351869631 * temp + 8.37935977E-5 * pow(temp,2.0) - 6.67125E-8 * pow(temp,3.0);

		}else if (temp < 273){
			dynamicViscosityTemp = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (temp > 433){
			dynamicViscosityTemp = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}


	}break;
	case ConvectionConfiguration::FLUID::TRANSFORMER_OIL:{

		if ((extTemp >=223) && (extTemp <= 373)){
			rho = 1055.04607 - 0.581753034 * extTemp - 6.40531689E-5 * pow(extTemp,2.0);
		}else if (extTemp < 223){
			rho =  1055.04607 - 0.581753034 * 223.0 - 6.40531689E-5 * pow(223.0,2.0);
		}else if (extTemp > 373){
			rho = 1055.04607 - 0.581753034 * 373.0 - 6.40531689E-5 * pow(373.0,2.0);
		}

		if ((extTemp >= 273) && (extTemp <= 433)){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * extTemp;
		}else if (extTemp < 273){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * 223.0;
		}else if (extTemp > 433){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * 373.0;
		}


		if ((extTemp >= 223) && (extTemp <= 293)){
			heatCapacity = -117056.38 + 1816.76208 * extTemp - 10.305786 * pow(extTemp,2.0) + 0.0256691919 * pow(extTemp,3.0) - 2.36742424E-5 * pow(extTemp,4.0);
		}else if ((extTemp > 293) && (extTemp <= 373 )){
			heatCapacity = -13408.1491 + 123.044152 * extTemp - 0.335401786 * pow(extTemp,2.0) + 3.125E-4 * pow(extTemp,3.0);
		}else if (extTemp < 223){
			heatCapacity = -117056.38 + 1816.76208 * 223.0 - 10.305786 * pow(223.0,2.0) + 0.0256691919 * pow(223.0,3.0) - 2.36742424E-5 * pow(223.0,4.0);
		}else if (extTemp > 373){
			heatCapacity = -13408.1491 + 123.044152 * 373.0 - 0.335401786 * pow(373.0,2.0) + 3.125E-4 * pow(373.0,3.0);
		}


		if ((extTemp >= 243) && (extTemp <= 273)){
			dynamicViscosity = 4492.20229 - 64.7408879 * extTemp + 0.349900959 * pow(extTemp,2.0) - 8.40477E-4 * pow(extTemp,3.0) + 7.57041667E-7 * pow(extTemp,4.0);
		}else if ((extTemp > 273) && (extTemp <= 373 )){
			dynamicViscosity = 91.4524999 - 1.33227058 * extTemp + 0.00777680216 * pow(extTemp,2.0) - 2.27271368E-5 *  pow(extTemp,3.0) + 3.32419673E-8 * pow(extTemp,4.0) - 1.94631023E-11 * pow(extTemp,5.0);
		}else if (extTemp < 243){
			dynamicViscosity = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (extTemp > 373){
			dynamicViscosity = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

		if ((temp >= 243) && (temp <= 273)){
			dynamicViscosityTemp = 4492.20229 - 64.7408879 * temp + 0.349900959 * pow(temp,2.0) - 8.40477E-4 * pow(temp,3.0) + 7.57041667E-7 * pow(temp,4.0);
		}else if ((temp > 273) && (temp <= 373 )){
			dynamicViscosityTemp = 91.4524999 - 1.33227058 * temp + 0.00777680216 * pow(temp,2.0) - 2.27271368E-5 *  pow(temp,3.0) + 3.32419673E-8 * pow(temp,4.0) - 1.94631023E-11 * pow(temp,5.0);
		}else if (temp < 243){
			dynamicViscosityTemp = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (temp > 373){
			dynamicViscosityTemp = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

	}break;
	default:
		eslog::error("Invalid convection fluid type.");
	}
}

double HeatTransferKernel::convectionHTC(
		const ConvectionConfiguration &convection,
		esint csize, double *coordinates, double time, double temp)
{
	double htc = 0;
	switch (convection.type) {
	case ConvectionConfiguration::TYPE::USER:
		convection.heat_transfer_coefficient.evaluator->evalVector(1, csize, coordinates, &temp, time, &htc);
		return htc;

	case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL: {
		double avgTemp, g, rho, dynamicViscosity, heatCapacity, thermalConductivity, dynamicViscosityTemp, extTemp;

		convection.external_temperature.evaluator->evalVector(1, csize, coordinates, &temp, time, &extTemp);
		avgTemp = (extTemp + temp) / 2;
		g = 9.81;

		convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);

		switch (convection.variant) {
		case ConvectionConfiguration::VARIANT::INCLINED_WALL: {
			double wallHeight, tiltAngle;
			convection.wall_height.evaluator->evalVector(1, csize, coordinates, &temp, time, &wallHeight);
			double RaL = pow(rho,2) * g * (1 / avgTemp) * heatCapacity * std::fabs(temp - extTemp  ) * pow(wallHeight, 3.0) / ( thermalConductivity * dynamicViscosity);
			convection.tilt_angle.evaluator->evalVector(1, csize, coordinates, &temp, time, &tiltAngle);
			tiltAngle *= M_PI / 180.0;
			if (RaL <= 1e9) {
				htc = (thermalConductivity / wallHeight) * (0.68 + (0.67 * cos(tiltAngle) * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermalConductivity / wallHeight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::VERTICAL_WALL: {
			double wallHeight;
			convection.wall_height.evaluator->evalVector(1, csize, coordinates, &temp, time, &wallHeight);
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(wallHeight, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (RaL <= 1e9) {
				htc = (thermalConductivity / wallHeight) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermalConductivity / wallHeight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_UP: {
			double length;
			convection.length.evaluator->evalVector(1, csize, coordinates, &temp, time, &length);
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(length, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (temp > extTemp) {
				if (RaL <= 1e7) {
					htc = thermalConductivity / length * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermalConductivity / length * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermalConductivity / length * 0.27 * pow(RaL,0.25);
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_DOWN: {
			double length;
			convection.length.evaluator->evalVector(1, csize, coordinates, &temp, time, &length);
			double RaL = pow(rho,2)	* g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) *pow(length, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (temp <= extTemp) {
				if (RaL <= 1e7) {
					htc = thermalConductivity / length * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermalConductivity / length * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermalConductivity / length * 0.27 * pow(RaL,0.25);
			}
		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_CYLINDER:{
			double diameter;
			convection.diameter.evaluator->evalVector(1, csize, coordinates, &temp, time, &diameter);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
			double Pr = dynamicViscosity * heatCapacity / thermalConductivity;

			if ( RaD > 10e12 ){
				// warning!!!!
				eslog::error("Validated only for RaD <= 10e12.");
			}
			htc = thermalConductivity / diameter * pow( 0.6 + ( 0.387*pow(RaD,1.0/6.0)/ pow( 1 + pow( 0.559/Pr, 9.0/16.0), 8.0/27.0) ) ,2.0);

		} break;
		case ConvectionConfiguration::VARIANT::SPHERE: {
			double diameter;
			convection.diameter.evaluator->evalVector(1, csize, coordinates, &temp, time, &diameter);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs( temp - extTemp) * pow(diameter, 3.0) / ( thermalConductivity * dynamicViscosity);
			double Pr = dynamicViscosity * heatCapacity / thermalConductivity;

			if ( RaD > 10e11 || Pr < 0.7 ){
				// warning!!!!
				eslog::error("Validated only for RaD <= 10e11 and Pr >= 0.7.");
			}
			htc = thermalConductivity / diameter * pow(2.0 + ( 0.589*pow(RaD,0.25)/ pow( 1 + pow( 0.469/Pr, 9.0/16.0), 4.0/9.0) ) ,2.0);
		} break;
		default:
			eslog::error("Invalid convection variant for EXTERNAL_NATURAL.");
		}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:{

		double avgTemp, g, rho, dynamicViscosity, heatCapacity, thermalConductivity,dynamicViscosityTemp, extTemp;

		convection.length.evaluator->evalVector(1, csize, coordinates, &temp, time, &extTemp);
		avgTemp = (extTemp + temp) / 2.0;
		g = 9.81;

		convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity );

		switch (convection.variant) {
		case ConvectionConfiguration::VARIANT::PARALLEL_PLATES: {
			double wallHeight, length;
			convection.wall_height.evaluator->evalVector(1, csize, coordinates, &temp, time, &wallHeight);
			convection.length.evaluator->evalVector(1, csize, coordinates, &temp, time, &length);
			double H_L = wallHeight / length;
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(length,3.0)/ ( thermalConductivity * dynamicViscosity);

			if (( RaL < H_L ) && (temp >  extTemp)) {
				htc = thermalConductivity / wallHeight * ( 1.0 / 24.0 ) * RaL;
			} else {
				double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs( temp - extTemp) * pow(length,3.0)/ ( thermalConductivity * dynamicViscosity);
				if (RaL <= 1e9) {
					htc = (thermalConductivity / length) * (0.68 + (0.67 * pow(RaL,0.25)) / (pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermalConductivity / length) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		case ConvectionConfiguration::VARIANT::CIRCULAR_TUBE: {
			double diameter, wallHeight;
			convection.diameter.evaluator->evalVector(1, csize, coordinates, &temp, time, &diameter);
			convection.wall_height.evaluator->evalVector(1, csize, coordinates, &temp, time, &wallHeight);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
			double H_D = wallHeight / diameter;

			if ( RaD < H_D ){
				htc = thermalConductivity / wallHeight * ( 1.0 / 128.0 ) * RaD;
			}else{

				double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
				if (RaD <= 1e9) {
					htc = (thermalConductivity / diameter) * (0.68 + (0.67 * pow(RaD,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermalConductivity / diameter) * pow(0.825 + (0.387 * pow(RaD,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		default:
			eslog::error("Invalid convection variant for INTERNAL_NATURAL.\n");
		}
	} break;

	case ConvectionConfiguration::TYPE::EXTERNAL_FORCED: {

			switch (convection.variant) {
			case ConvectionConfiguration::VARIANT::AVERAGE_PLATE: {

				double avgTemp, rho, dynamicViscosity, heatCapacity, thermalConductivity,dynamicViscosityTemp, extTemp, length, fluidVelocity;
				convection.external_temperature.evaluator->evalVector(1, csize, coordinates, &temp, time, &extTemp);
				convection.length.evaluator->evalVector(1, csize, coordinates, &temp, time, &length);
				convection.fluid_velocity.evaluator->evalVector(1, csize, coordinates, &temp, time, &fluidVelocity);

				avgTemp = (extTemp + temp) / 2.0;

				convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);


				double Re = rho * fluidVelocity * length / dynamicViscosity;
				double Pr = dynamicViscosity * heatCapacity / thermalConductivity;
				if (Re <= 5e5) {
					htc = 2 * (thermalConductivity / length) * ((0.3387 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.5)) / (pow(1 + pow(0.0468 / Pr, 2.0 / 3.0), 0.25)));
				} else {
					htc = 2 * (thermalConductivity / length) * pow(Pr, 1.0 / 3.0)	* (0.037 * pow(Re, 0.8) - 871);
				}

			} break;
			default:
				eslog::error("Invalid convection variant for EXTERNAL_FORCED.\n");
			}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_FORCED:{

			switch (convection.variant) {
			case ConvectionConfiguration::VARIANT::TUBE: {

				double extTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity, fluidVelocity, diameter;

				convection.external_temperature.evaluator->evalVector(1, csize, coordinates, &temp, time, &extTemp);
				convection.fluid_velocity.evaluator->evalVector(1, csize, coordinates, &temp, time, &fluidVelocity);
				convection.diameter.evaluator->evalVector(1, csize, coordinates, &temp, time, &diameter);

				convectionMaterialParameters(convection, csize, coordinates, time, temp, extTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);

				double Re = rho * fluidVelocity * diameter / dynamicViscosity;
				double Pr = dynamicViscosity * heatCapacity / thermalConductivity;
				double n = temp < extTemp ? 0.3 : 0.4;
				htc = thermalConductivity / diameter;
				if (Re <= 2500) {
					htc *= 3.66;
				} else {
					htc *= 0.027 * pow(Re, .8) * pow(Pr, n)	* pow(dynamicViscosity / dynamicViscosityTemp, 0.14);
				}
			}break;
			default:
				eslog::error("Invalid convection variant for EXTERNAL_FORCED.\n");
			}
	}break;

	default:
		eslog::error("Invalid convection TYPE.\n");
	}

	return htc;
}
