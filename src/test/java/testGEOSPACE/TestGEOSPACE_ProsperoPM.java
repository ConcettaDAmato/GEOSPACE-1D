/*
* GNU GPL v3 License
 *
 * Copyright 2019 Concetta D'Amato
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package testGEOSPACE;
import java.net.URISyntaxException;
import java.util.*;
import org.junit.Test;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.hortonmachine.gears.io.rasterreader.OmsRasterReader;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
//import org.jgrasstools.gears.libs.monitor.PrintStreamProgressMonitor;

import it.geoframe.blogspot.brokergeo.solver.*;
import it.geoframe.blogspot.brokergeo.data.*;
import it.geoframe.blogspot.buffer.buffertowriter.*;
import it.geoframe.blogspot.whetgeo1d.richardssolver.*;
import it.geoframe.blogspot.netcdf.monodimensionalproblemtimedependent.*;
import it.geoframe.blogspot.geoet.inout.InputReaderMain;
import it.geoframe.blogspot.geoet.inout.OutputWriterMain;
import it.geoframe.blogspot.geoet.transpiration.solver.ProsperoSolverMain;
import it.geoframe.blogspot.geoet.rootdensity.solver.RootDensitySolverMain;
import it.geoframe.blogspot.geoet.soilevaporation.solver.PMEvaporationFromSoilCanopySolverMain;
import it.geoframe.blogspot.geoet.stressfactor.solver.*;
import it.geoframe.blogspot.geoet.totalEvapoTranspiration.TotalEvapoTranspirationSolverMain;


/**
 * Test  GEOSPACE
 * This is the test for the GEOSPACE MODEL
 * @author Concetta D'Amato, Niccolo' Tubini, Michele Botazzi and Riccardo Rigon.  
 */
public class TestGEOSPACE_ProsperoPM {

	@Test
	public void Test() throws Exception {
		
		String startDate= "2015-04-01 00:00";
        String endDate	= "2015-06-01 02:00";
        String fId = "ID";
        String Id = "1";
        String site = "Cavone/";
		int timeStepMinutes = 60;
		String lab = "waterstress"; ////richards - potential - waterstress -  environmentalstress - totalstress - potential_evaporation
		String lab2 = "test02";
		
		
     	String pathTopBC    ="data/"+site+Id+"/precip_1.csv";
		String pathBottomBC ="data/"+site+Id+"/Cavone_0.csv";
		String pathGrid     ="data/Grid_NetCDF/Grid_GEOSPACE_test.nc";
		String pathSaveDates="data/"+site+Id+"/Cavone_1.csv"; 
		String pathOutput = "output/GEOSPACE/JavaTestGEOSPACE_"+lab+lab2+".nc";
		String outputDescription = "\n"
				+ "Initial condition hydrostatic no ponding\n		"
				+ "Bottom Dirichlet\n		"
				+ "Grid input file: " + pathGrid +"\n		"
				+ "TopBC input file: " + pathTopBC +"\n		"
				+ "BottomBC input file: " + pathBottomBC +"\n		"
				+ "DeltaT: 50s\n		"
				+ "Picard iteration: 1\n		"
			    + "Interface k: max";
		
		String topBC = "Top Coupled";
		String bottomBC = "Bottom Free drainage";
		
       
        OmsRasterReader DEMreader = new OmsRasterReader();
		DEMreader.file = "data/"+site+Id+"/dem_1.tif";
		//DEMreader.fileNovalue = -9999.0;
		//DEMreader.geodataNovalue = Double.NaN;
		DEMreader.process();
		GridCoverage2D digitalElevationModel = DEMreader.outRaster;
		
		String inPathToTemperature 				="data/"+site+Id+"/airT_1.csv";
        String inPathToWind 					="data/"+site+Id+"/Wind_1.csv";
        String inPathToRelativeHumidity 		="data/"+site+Id+"/RH_1.csv";
        String inPathToShortWaveRadiationDirect ="data/"+site+Id+"/ShortwaveDirect_1.csv";
        String inPathToShortWaveRadiationDiffuse="data/"+site+Id+"/ShortwaveDiffuse_1.csv";
        String inPathToLWRad 					="data/"+site+Id+"/LongDownwelling_1.csv";
        String inPathToNetRad 					="data/"+site+Id+"/Net_1.csv";
        String inPathToSoilHeatFlux 			="data/"+site+Id+"/GHF_all_1.csv";
        String inPathToPressure 				="data/"+site+Id+"/Pres_1.csv";
        String inPathToLai 						="data/"+site+Id+"/LAI_4.csv";
        String inPathToCentroids 				="data/"+site+Id+"/centroids_ID_1.shp";
        
        //String inPathToRoot 					="data/"+site+Id+"/RootDepth_2306_04.csv";
        //String inPathToCanopyHeight				="data/"+site+Id+"/CanopyH_2206_09.csv";
       
        
        
        String outPathToLatentHeatSun			="output/GEOSPACE/LatentHeatSun_"+lab+"_"+lab2+".csv";
        String outPathToLatentHeatShadow		="output/GEOSPACE/LatentHeatShadow_"+lab+"_"+lab2+".csv";
        String outPathToSoilFluxEvaporation		="output/GEOSPACE/FluxEvaporation_"+lab+"_"+lab2+".csv";
        String outPathToFluxTranspiration		="output/GEOSPACE/FluxTranspiration_"+lab+"_"+lab2+".csv";
        String outPathToFluxEvapoTranspiration	="output/GEOSPACE/FluxEvapoTranspiration_"+lab+"_"+lab2+".csv";
        String outPathToTranspiration			="output/GEOSPACE/Transpiration_"+lab+"_"+lab2+".csv";
        String outPathToEvapoTranspiration		="output/GEOSPACE/EvapoTranspiration_"+lab+"_"+lab2+".csv";
        String outPathToSoilEvaporation 		="output/GEOSPACE/Evaporation_"+lab+"_"+lab2+".csv";
        
		String outPathToLeafTemperatureSun		="output/GEOSPACE/LeafTemperatureSun_"+lab+"_"+lab2+".csv";
		String outPathToLeafTemperatureShadow	="output/GEOSPACE/LeafTemperatureSh_"+lab+"_"+lab2+".csv";
		String outPathToSensibleSun				="output/GEOSPACE/sensibleSun_"+lab+"_"+lab2+".csv";
		String outPathToSensibleShadow			="output/GEOSPACE/sensibleShadow_"+lab+"_"+lab2+".csv";
        String outPathToRadiationSoil 			="output/GEOSPACE/RadiationSoil_"+lab+"_"+lab2+".csv";
		String outPathToRadiationSun			="output/GEOSPACE/RadSun_"+lab+"_"+lab2+".csv";
		String outPathToRadiationShadow			="output/GEOSPACE/RadShadow_"+lab+"_"+lab2+".csv";
		String outPathToCanopy					="output/GEOSPACE/Canopy_"+lab+"_"+lab2+".csv";
		String outPathToVPD						="output/GEOSPACE/VPD_"+lab+"_"+lab2+".csv";
		
		RichardsRootSolver1DMain R1DSolver     		= new RichardsRootSolver1DMain();
		GEOSPACEBuffer1D buffer 		  			= new GEOSPACEBuffer1D();
		WriteNetCDFGEOSPACE1DDouble writeNetCDF 	= new WriteNetCDFGEOSPACE1DDouble();
		ReadNetCDFGEOSPACEGrid1D readNetCDF    		= new ReadNetCDFGEOSPACEGrid1D();	
		
		InputReaderMain Input 											= new InputReaderMain();
		OutputWriterMain Output 										= new OutputWriterMain();
		
		JarvisStressFactorSolverMain JarvisStressFactor   				= new JarvisStressFactorSolverMain();
		ProsperoSolverMain Prospero 									= new ProsperoSolverMain();
		PMEvaporationFromSoilCanopySolverMain PMsoilevaporation 		= new PMEvaporationFromSoilCanopySolverMain();
		TotalEvapoTranspirationSolverMain TotalEvapoTranspiration 		= new TotalEvapoTranspirationSolverMain();
		
		InputDataMain InputBroker										= new InputDataMain();
		ETsBrokerTwoFluxesSolverMain ETsBrokerSolver 					= new ETsBrokerTwoFluxesSolverMain(); 
		RootDensitySolverMain RootDensitySolver							= new RootDensitySolverMain();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// model's variables /////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int writeFrequency = 1;
		
		R1DSolver.stationID = 1;
		R1DSolver.beta0 = -766.45;
		R1DSolver.referenceTemperatureSWRC = 293.15;
		R1DSolver.maxPonding = 0;
		
		R1DSolver.typeClosureEquation = new String[] {"Water Depth", "Van Genuchten"};
		R1DSolver.typeEquationState = new String[] {"Water Depth", "Van Genuchten"};
		
		R1DSolver.typeUHCModel = new String[] {"", "Mualem Van Genuchten"};
		R1DSolver.typeUHCTemperatureModel = "notemperature"; //"Ronan1998";
		R1DSolver.interfaceHydraulicConductivityModel = "max";
		R1DSolver.delta = 0;
		R1DSolver.tTimeStep = 3600;
		R1DSolver.timeDelta = 3600;
		R1DSolver.newtonTolerance = Math.pow(10,-9);
		R1DSolver.nestedNewton =1;
		R1DSolver.picardIteration = 1;

		writeNetCDF.outVariables = new String[] {"darcyVelocity"};
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.soilHydraulicConductivityModel = "Mualem VG no temperature";
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.fileSizeMax = 10000;
		
		Input.idCentroids="ID";
		Input.centroidElevation="Elevation";
		
		
		//----------------------Canopy----------------------
		Input.canopyHeight = 2.5;
		//Input.canopyHeightType = "canopyHeightGrowth";
		Prospero.typeOfCanopy="multilayer";
		
		//----------------------Root----------------------
		Input.rootDepth  = -2;
		//Input.rootType = "rootGrowth"; //variabile che controlla la lettura di una serie temporale per la profondità delle radici
		//RootDensitySolver.rootDensityModel = "ExponentialGrowthMethod"; //CostantMethod, CostantGrowthMethod, LinearGrowthMethod, ExponentialGrowthMethod ----------------//metodi per il calcolo della densità ad ogni volume di controllo
		//Input.growthRateRoot= 0.001;
				
		//----------------------Evaporation Layer----------------------
		JarvisStressFactor.etaE  = -0.2; //depth of the evaporation layer
		
		
		//----------------------Stress Factor----------------------
		JarvisStressFactor.stressFactorModel = "LinearStressFactor";
		JarvisStressFactor.representativeTranspirationSFModel = "RootDensityWeightedMethod"; //SizeWightedMethod, AverageMethod, RootDensityWeightedMethod
		JarvisStressFactor.representativeEvaporationSFModel = "AverageMethod"; //SizeWightedMethod, AverageMethod, RootDensityWeightedMethod
		
		//----------------------EvapoTranspiration Splitter Methods----------------------
		ETsBrokerSolver.representativeEsModel = "AverageWaterWeightedMethod";  	//SizeWaterWeightedMethod, AverageWaterWeightedMethod //SizeWightedMethod, AverageWeightedMethod
		ETsBrokerSolver.representativeTsModel = "RootWaterWeightedMethod"; //SizeWaterWeightedMethod, AverageWaterWeightedMethod, RootWaterWeightedMethod //SizeWightedMethod, AverageWeightedMethod, RootWeightedMethod
       	
		
		
		
		JarvisStressFactor.useRadiationStress    = false;
		JarvisStressFactor.useTemperatureStress  = false;
		JarvisStressFactor.useVDPStress		   = false;	
		JarvisStressFactor.useWaterStress = true;
		//JarvisStressFactor.defaultStress      = 1.0;
						
		JarvisStressFactor.alpha = 0.005;
		JarvisStressFactor.thetaR = 0.9;
		JarvisStressFactor.VPD0 = 5.0;
		JarvisStressFactor.Tl = -5.0;
		JarvisStressFactor.T0 = 20.0;
		JarvisStressFactor.Th = 45.0;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
		OmsTimeSeriesIteratorReader topBCReader    = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		OmsShapefileFeatureReader centroidsReader = new OmsShapefileFeatureReader();
        centroidsReader.file = inPathToCentroids;
		centroidsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = centroidsReader.geodata;
		Input.inCentroids = stationsFC;
		Input.inDem = digitalElevationModel;	
        OmsTimeSeriesIteratorReader temperatureReader		= getTimeseriesReader(inPathToTemperature, fId, startDate, endDate, timeStepMinutes);
        OmsTimeSeriesIteratorReader windReader 		 		= getTimeseriesReader(inPathToWind, fId, startDate, endDate, timeStepMinutes);
        OmsTimeSeriesIteratorReader humidityReader 			= getTimeseriesReader(inPathToRelativeHumidity, fId, startDate, endDate, timeStepMinutes);
        OmsTimeSeriesIteratorReader shortwaveReaderDirect 	= getTimeseriesReader(inPathToShortWaveRadiationDirect, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader shortwaveReaderDiffuse 	= getTimeseriesReader(inPathToShortWaveRadiationDiffuse, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader longwaveReader 			= getTimeseriesReader(inPathToLWRad, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader pressureReader 			= getTimeseriesReader(inPathToPressure, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader leafAreaIndexReader		= getTimeseriesReader(inPathToLai, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader soilHeatFluxReader 		= getTimeseriesReader(inPathToSoilHeatFlux, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader netRadReader 			= getTimeseriesReader(inPathToNetRad, fId, startDate, endDate,timeStepMinutes);
		
        //OmsTimeSeriesIteratorReader rootReader 			= getTimeseriesReader(inPathToRoot, fId, startDate, endDate,timeStepMinutes);
        //OmsTimeSeriesIteratorReader canopyHeightReader 			= getTimeseriesReader(inPathToCanopyHeight, fId, startDate, endDate,timeStepMinutes);
        
        
        OmsTimeSeriesIteratorWriter latentHeatSunWriter = new OmsTimeSeriesIteratorWriter();
		latentHeatSunWriter.file = outPathToLatentHeatSun;
		latentHeatSunWriter.tStart = startDate;
		latentHeatSunWriter.tTimestep = timeStepMinutes;
		latentHeatSunWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter latentHeatShadowWriter = new OmsTimeSeriesIteratorWriter();
		latentHeatShadowWriter.file = outPathToLatentHeatShadow;
		latentHeatShadowWriter.tStart = startDate;
		latentHeatShadowWriter.tTimestep = timeStepMinutes;
		latentHeatShadowWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter FluxTranspirationWriter = new OmsTimeSeriesIteratorWriter();
		FluxTranspirationWriter.file = outPathToFluxTranspiration;
		FluxTranspirationWriter.tStart = startDate;
		FluxTranspirationWriter.tTimestep = timeStepMinutes;
		FluxTranspirationWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter FluxEvaporationWriter = new OmsTimeSeriesIteratorWriter();
		FluxEvaporationWriter.file = outPathToSoilFluxEvaporation;
		FluxEvaporationWriter.tStart = startDate;
		FluxEvaporationWriter.tTimestep = timeStepMinutes;
		FluxEvaporationWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter FluxEvapoTranspirationWriter = new OmsTimeSeriesIteratorWriter();
		FluxEvapoTranspirationWriter.file = outPathToFluxEvapoTranspiration;
		FluxEvapoTranspirationWriter.tStart = startDate;
		FluxEvapoTranspirationWriter.tTimestep = timeStepMinutes;
		FluxEvapoTranspirationWriter.fileNovalue="-9999";

		OmsTimeSeriesIteratorWriter EvapoTranspirationWriter = new OmsTimeSeriesIteratorWriter();
		EvapoTranspirationWriter.file = outPathToEvapoTranspiration;
		EvapoTranspirationWriter.tStart = startDate;
		EvapoTranspirationWriter.tTimestep = timeStepMinutes;
		EvapoTranspirationWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter TranspirationWriter = new OmsTimeSeriesIteratorWriter();
		TranspirationWriter.file = outPathToTranspiration;
		TranspirationWriter.tStart = startDate;
		TranspirationWriter.tTimestep = timeStepMinutes;
		TranspirationWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter EvaporationWriter = new OmsTimeSeriesIteratorWriter();
		EvaporationWriter.file = outPathToSoilEvaporation;
		EvaporationWriter.tStart = startDate;
		EvaporationWriter.tTimestep = timeStepMinutes;
		EvaporationWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter leafTemperatureSunWriter = new OmsTimeSeriesIteratorWriter();
		leafTemperatureSunWriter.file = outPathToLeafTemperatureSun;
		leafTemperatureSunWriter.tStart = startDate;
		leafTemperatureSunWriter.tTimestep = timeStepMinutes;
		leafTemperatureSunWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter leafTemperatureShadowWriter = new OmsTimeSeriesIteratorWriter();
		leafTemperatureShadowWriter.file = outPathToLeafTemperatureShadow;
		leafTemperatureShadowWriter.tStart = startDate;
		leafTemperatureShadowWriter.tTimestep = timeStepMinutes;
		leafTemperatureShadowWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter radiationSunWriter = new OmsTimeSeriesIteratorWriter();
		radiationSunWriter.file = outPathToRadiationSun;
		radiationSunWriter.tStart = startDate;
		radiationSunWriter.tTimestep = timeStepMinutes;
		radiationSunWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter radiationShadowWriter = new OmsTimeSeriesIteratorWriter();
		radiationShadowWriter.file = outPathToRadiationShadow;
		radiationShadowWriter.tStart = startDate;
		radiationShadowWriter.tTimestep = timeStepMinutes;
		radiationShadowWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter sensibleSunWriter = new OmsTimeSeriesIteratorWriter();
		sensibleSunWriter.file = outPathToSensibleSun;
		sensibleSunWriter.tStart = startDate;
		sensibleSunWriter.tTimestep = timeStepMinutes;
		sensibleSunWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter sensibleShadowWriter = new OmsTimeSeriesIteratorWriter();
		sensibleShadowWriter.file = outPathToSensibleShadow;
		sensibleShadowWriter.tStart = startDate;
		sensibleShadowWriter.tTimestep = timeStepMinutes;
		sensibleShadowWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter radiationSoilWriter = new OmsTimeSeriesIteratorWriter();
		radiationSoilWriter.file = outPathToRadiationSoil;
		radiationSoilWriter.tStart = startDate;
		radiationSoilWriter.tTimestep = timeStepMinutes;
		radiationSoilWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter canopyWriter = new OmsTimeSeriesIteratorWriter();
		canopyWriter.file = outPathToCanopy;
		canopyWriter.tStart = startDate;
		canopyWriter.tTimestep = timeStepMinutes;
		canopyWriter.fileNovalue="-9999";
		
		OmsTimeSeriesIteratorWriter vapourPressureDeficitWriter = new OmsTimeSeriesIteratorWriter();
		vapourPressureDeficitWriter.file = outPathToVPD;
		vapourPressureDeficitWriter.tStart = startDate;
		vapourPressureDeficitWriter.tTimestep = timeStepMinutes;
		vapourPressureDeficitWriter.fileNovalue="-9999";
		
		readNetCDF.richardsGridFilename = pathGrid;
		readNetCDF.read();
		
		R1DSolver.z = readNetCDF.z;
		R1DSolver.spaceDeltaZ = readNetCDF.spaceDelta;
		R1DSolver.psiIC = readNetCDF.psiIC;
		R1DSolver.temperature = readNetCDF.temperature;
		R1DSolver.controlVolume = readNetCDF.controlVolume;
		R1DSolver.ks = readNetCDF.Ks;
		R1DSolver.thetaS = readNetCDF.thetaS;
		R1DSolver.thetaR = readNetCDF.thetaR;
		R1DSolver.par1SWRC = readNetCDF.par1SWRC;
		R1DSolver.par2SWRC = readNetCDF.par2SWRC;
		R1DSolver.par3SWRC = readNetCDF.par3SWRC;
		R1DSolver.par4SWRC = readNetCDF.par4SWRC;
		R1DSolver.par5SWRC = readNetCDF.par5SWRC;
		R1DSolver.alphaSpecificStorage = readNetCDF.alphaSS;
		R1DSolver.betaSpecificStorage = readNetCDF.betaSS;
		R1DSolver.inEquationStateID = readNetCDF.equationStateID;
		R1DSolver.inParameterID = readNetCDF.parameterID;
		R1DSolver.thetaWP = readNetCDF.thetaWP;
		R1DSolver.thetaFC = readNetCDF.thetaFC;
		R1DSolver.topBCType = topBC;
		R1DSolver.bottomBCType = bottomBC;
		
		buffer.writeFrequency = writeFrequency;
	
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathRichardsBottomBC = pathBottomBC; 
		writeNetCDF.pathRichardsTopBC = pathTopBC; 
		writeNetCDF.bottomRichardsBC = bottomBC;
		writeNetCDF.topRichardsBC = topBC;
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psiIC = readNetCDF.psiIC;
		writeNetCDF.temperature = readNetCDF.temperature;
		writeNetCDF.rootIC = readNetCDF.rootIC;

		writeNetCDF.timeUnits = "Minutes since 01/01/1970 01:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		
				
		JarvisStressFactor.thetaWp	= readNetCDF.thetaWP;
		JarvisStressFactor.thetaFc 	= readNetCDF.thetaFC;
		JarvisStressFactor.ID      	= readNetCDF.parameterID;
		JarvisStressFactor.deltaZ  	= readNetCDF.spaceDelta;
		InputBroker.deltaZ 			= readNetCDF.spaceDelta;
		JarvisStressFactor.z      	= readNetCDF.z;
		InputBroker.z      			= readNetCDF.z;
		
		//InputBroker.etaR   			= JarvisStressFactor.etaR;
		InputBroker.etaE   			= JarvisStressFactor.etaE;
		//InputBroker.rootIC          = readNetCDF.rootIC;
		Input.z      			= readNetCDF.z;
		Input.rootIC          = readNetCDF.rootIC;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
//////////////////////// START /////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
			while( topBCReader.doProcess) {
	        
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			R1DSolver.inTopBC= bCValueMap;
			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			R1DSolver.inBottomBC = bCValueMap;
			
			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			R1DSolver.inSaveDate = bCValueMap;
			
			R1DSolver.inCurrentDate = topBCReader.tCurrent;

	 
			temperatureReader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = temperatureReader.outData;
			Input.inAirTemperature = id2ValueMap;
			//Input.doHourly = true;
            Output.doFullPrint = false;
            Input.tStartDate = startDate;
            Input.temporalStep = timeStepMinutes;

            windReader.nextRecord();
            id2ValueMap = windReader.outData;
            Input.inWindVelocity = id2ValueMap;
            
            humidityReader.nextRecord();
            id2ValueMap = humidityReader.outData;
            Input.inRelativeHumidity = id2ValueMap;
            
            shortwaveReaderDirect.nextRecord();
            id2ValueMap = shortwaveReaderDirect.outData;
            Input.inShortWaveRadiationDirect = id2ValueMap; 
            
            shortwaveReaderDiffuse.nextRecord();
            id2ValueMap = shortwaveReaderDiffuse.outData;
            Input.inShortWaveRadiationDiffuse = id2ValueMap;
            
            longwaveReader.nextRecord();
            id2ValueMap = longwaveReader.outData;
            Input.inLongWaveRadiation = id2ValueMap;
            
            soilHeatFluxReader.nextRecord();
            id2ValueMap = soilHeatFluxReader.outData;
            Input.inSoilFlux = id2ValueMap;
            
            pressureReader.nextRecord();
            id2ValueMap = pressureReader.outData;
            Input.inAtmosphericPressure = id2ValueMap;
            
            leafAreaIndexReader.nextRecord();
            id2ValueMap = leafAreaIndexReader.outData;
            Input.inLeafAreaIndex = id2ValueMap;
            
            netRadReader.nextRecord();
            id2ValueMap = netRadReader.outData;
            Input.inNetRadiation = id2ValueMap;
            
            //rootReader.nextRecord();
            //id2ValueMap = rootReader.outData;
            //Input.inRootDepth = id2ValueMap;
            
            //canopyHeightReader.nextRecord();
            //id2ValueMap = canopyHeightReader.outData;
            //Input.inCanopyHeight = id2ValueMap;

/////////////////////////////// SOLVE /////////////////////////
			
            R1DSolver.solve();

			
			
			JarvisStressFactor.theta = R1DSolver.thetasNew;
			
			
			Input.process();
			InputBroker.etaR = Input.defRootDepth;
			
			RootDensitySolver.solve();
			
			InputBroker.rootDensity = RootDensitySolver.defRootDensity;
		
			
			JarvisStressFactor.solve();
            
    
			PMsoilevaporation.evaporationStressWater = JarvisStressFactor.evaporationStressWater;
            Prospero.stressSun = JarvisStressFactor.stressSun;
            Prospero.stressShade = JarvisStressFactor.stressShade;
            
            Prospero.process();
            
            PMsoilevaporation.process();
            
            TotalEvapoTranspiration.transpiration = Prospero.transpiration;
            TotalEvapoTranspiration.evaporation = PMsoilevaporation.evaporation;
            TotalEvapoTranspiration.process();
            
            InputBroker.g = JarvisStressFactor.g;
            InputBroker.GnT = JarvisStressFactor.GnT;
            InputBroker.GnE = JarvisStressFactor.GnE;
            InputBroker.transpiration = Prospero.transpiration;
            InputBroker.evaporation = PMsoilevaporation.evaporation;
           
            InputBroker.process();
			
            /////////////////////SOLO RICHARDS///////////////////////////////
			//InputBroker.transpiration = 0;
			//InputBroker.evaporation = 0;
			
            ETsBrokerSolver.useWaterStress=JarvisStressFactor.useWaterStress;
			ETsBrokerSolver.solve();
			
			R1DSolver.stressedETs = ETsBrokerSolver.StressedETs;
			
			buffer.inputDate = R1DSolver.inCurrentDate;
			buffer.doProcessBuffer = R1DSolver.doProcessBuffer;
			buffer.inputVariableRichards = R1DSolver.outputToBuffer;
			buffer.inputVariableBroker = ETsBrokerSolver.outputToBuffer;
			buffer.inputVariableStressFactor = JarvisStressFactor.outputToBuffer;
			
			
			buffer.solve();
			writeNetCDF.variables = buffer.myVariable;			
			writeNetCDF.doProcess = topBCReader.doProcess;
			
			writeNetCDF.writeNetCDF();
			
            Output.process();
			
			
            latentHeatSunWriter.inData = Output.outLatentHeatSun;
            latentHeatSunWriter.writeNextLine();

			latentHeatShadowWriter.inData = Output.outLatentHeatShade;
            latentHeatShadowWriter.writeNextLine();		
            
            FluxTranspirationWriter.inData = Output.outFluxTranspiration;
            FluxTranspirationWriter.writeNextLine();
            
            FluxEvaporationWriter.inData = Output.outFluxEvaporation;
            FluxEvaporationWriter.writeNextLine();
            
            FluxEvapoTranspirationWriter.inData = Output.outFluxEvapoTranspiration;
            FluxEvapoTranspirationWriter.writeNextLine();	
			
			EvapoTranspirationWriter.inData = Output.outEvapoTranspiration;
			EvapoTranspirationWriter.writeNextLine();
			
			TranspirationWriter.inData = Output.outTranspiration;
			TranspirationWriter.writeNextLine();
    		
			EvaporationWriter.inData = Output.outEvaporation;
			EvaporationWriter.writeNextLine();
			
			leafTemperatureSunWriter.inData = Output.outLeafTemperature;
			leafTemperatureSunWriter.writeNextLine();			 	

			leafTemperatureShadowWriter.inData = Output.outLeafTemperatureShade;
			leafTemperatureShadowWriter.writeNextLine();

			
			if (Output.doFullPrint == true) {
						 	

			radiationSunWriter.inData = Output.outRadiation;
			radiationSunWriter.writeNextLine();			 	

			radiationShadowWriter.inData = Output.outRadiationShade;
			radiationShadowWriter.writeNextLine();			 	
			
			sensibleSunWriter.inData = Output.outSensibleHeat;
			sensibleSunWriter.writeNextLine();			 	
			
			sensibleShadowWriter.inData = Output.outSensibleHeatShade;
			sensibleShadowWriter.writeNextLine();
			
			radiationSoilWriter.inData = Output.outRadiationSoil;
			radiationSoilWriter.writeNextLine();
			
			canopyWriter.inData = Output.outCanopy;
			canopyWriter.writeNextLine();
			
			vapourPressureDeficitWriter.inData = Output.outVapourPressureDeficit;
			vapourPressureDeficitWriter.writeNextLine();
			
			
			radiationSunWriter.close();
			radiationShadowWriter.close();
			sensibleSunWriter.close();
			sensibleShadowWriter.close();
			radiationSoilWriter.close();
			canopyWriter.close();
			vapourPressureDeficitWriter.close();
			}	
			
		
			
		}	
	   		
		topBCReader.close();
		bottomBCReader.close();

		temperatureReader.close();        
        windReader.close();
        humidityReader.close();     
        shortwaveReaderDirect.close();
        shortwaveReaderDiffuse.close();
        longwaveReader.close();
        soilHeatFluxReader.close();
        pressureReader.close();
        leafAreaIndexReader.close();
        
        //rootReader.close();
        //canopyHeightReader.close();
                
        latentHeatSunWriter.close();
		latentHeatShadowWriter.close();
		FluxEvapoTranspirationWriter.close();
		FluxEvaporationWriter.close();
		FluxTranspirationWriter.close();
		EvapoTranspirationWriter.close();
		TranspirationWriter.close();
		EvaporationWriter.close();
		leafTemperatureSunWriter.close();
		leafTemperatureShadowWriter.close();

	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
