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

import it.geoframe.blogspot.brokergeo.solver.*;
import it.geoframe.blogspot.brokergeo.data.*;
import it.geoframe.blogspot.buffer.buffertowriter.*;
import it.geoframe.blogspot.whetgeo1d.richardssolver.*;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.*;
import it.geoframe.blogspot.geoet.inout.InputReaderMain;
import it.geoframe.blogspot.geoet.inout.OutputWriterMain;
import it.geoframe.blogspot.geoet.priestleytaylor.PriestleyTaylorActualETSolverMain;
import it.geoframe.blogspot.geoet.rootdensity.solver.RootDensitySolverMain;
//import it.geoframe.blogspot.geoet.prospero.solver.ProsperoSolverMain;
//import it.geoframe.blogspot.geoet.soilevaporation.solver.PMEvaporationFromSoilAfterCanopySolverMain;
import it.geoframe.blogspot.geoet.stressfactor.solver.*;
//import it.geoframe.blogspot.geoet.totalEvapoTranspiration.TotalEvapoTranspirationSolverMain;


/**
 * Test  GEOSPCE
 * This is the test for the GEOSPACE MODEL
 * @author Concetta D'Amato, Niccolo' Tubini, Michele Botazzi and Riccardo Rigon.  
 */
public class TestGEOSPACE_PriestleyTaylor {

	@Test
	public void Test() throws Exception {
		
		String startDate= "2018-05-10 00:00";
        String endDate	= "2018-06-29 02:00";
        String fId = "ID";
        String Id = "1";
        String site = "SpikeII/";
		int timeStepMinutes = 60;
		String lab = "SpikePT"; ////richards - potential - waterstress -  environmentalstress - totalstress - potential_evaporation
		String lab2 = "testconsole";
		
		
     	String pathTopBC    ="data/"+site+"/Prec_Irrig_Height_hourly.csv";
		String pathBottomBC ="data/"+site+"/SpikeII_0.csv";
		String pathGrid     ="data/Grid_NetCDF/Grid_GEOSPACE_test.nc";
		String pathSaveDates="data/"+site+"/saveSpikeII_hourly.csv"; 
		String pathOutput = "output/GEOSPACE/Java_"+lab+"_"+lab2+".nc";
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
		String bottomBC = "Bottom Free Drainage";
		
        //PrintStreamProgressMonitor pm = new PrintStreamProgressMonitor(System.out, System.out);
        OmsRasterReader DEMreader = new OmsRasterReader();
		DEMreader.file = "data/"+site+"/DemSpikeIIcut.tif";
		//DEMreader.fileNovalue = -9999.0;
		//DEMreader.geodataNovalue = Double.NaN;
		DEMreader.process();
		GridCoverage2D digitalElevationModel = DEMreader.outRaster;
		
		String inPathToTemperature 					="data/"+site+"/AirTemperature_hourly.csv";
        //String inPathToWind 						="data/"+site+Id+"/Wind_1.csv";
        //String inPathToRelativeHumidity 			="data/"+site+Id+"/RH_1.csv";
		//String inPathToShortWaveRadiationDirect 	="data/"+site+Id+"/Cavone_ShortwaveDirect_1.csv";
		//String inPathToShortWaveRadiationDiffuse	="data/"+site+Id+"/Cavone_ShortwaveDiffuse_1.csv";
		//String inPathToLWRad 						="data/"+site+Id+"/Cavone_LongDownwelling_1.csv";
        String inPathToNetRad 						="data/"+site+"/Solar_Radiation_mean_hourly.csv";
        String inPathToSoilHeatFlux 				="data/"+site+"/SpikeII_nan.csv";
        String inPathToPressure 					="data/"+site+"/SpikeII_nan.csv";
        // String inPathToLai 						="data/"+site+Id+"/LAI_sin.csv";
        String inPathToCentroids 					="data/"+site+"/centroid.shp";
        
        //String outPathToLatentHeatSun				="output/GEOSPACE/LatentHeatSun_"+lab+"_"+lab2+".csv";
        //String outPathToLatentHeatShadow			="output/GEOSPACE/LatentHeatShadow_"+lab+"_"+lab2+".csv";
        //String outPathToSoilFluxEvaporation		="output/GEOSPACE/FluxEvaporation_"+lab+"_"+lab2+".csv";
        //String outPathToFluxTranspiration			="output/GEOSPACE/FluxTranspiration_"+lab+"_"+lab2+".csv";
        String outPathToFluxEvapoTranspiration		="output/GEOSPACE/FluxEvapoTranspiration_"+lab+"_"+lab2+".csv";
        //String outPathToTranspiration				="output/GEOSPACE/Transpiration_"+lab+"_"+lab2+".csv";
        String outPathToEvapoTranspiration			="output/GEOSPACE/EvapoTranspiration_"+lab+"_"+lab2+".csv";
        //String outPathToSoilEvaporation 			="output/GEOSPACE/Evaporation_"+lab+"_"+lab2+".csv";
        
        //String outPathToLeafTemperatureSun		="output/GEOSPACE/LeafTemperatureSun.csv";
        //String outPathToLeafTemperatureShadow		="output/GEOSPACE/LeafTemperatureSh.csv";
        //String outPathToSensibleSun				="output/GEOSPACE/sensibleSun.csv";
        //String outPathToSensibleShadow			="output/GEOSPACE/sensibleShadow.csv";
        //String outPathToRadiationSoil 			="output/GEOSPACE/RadiationSoil.csv";
        //String outPathToRadiationSun				="output/GEOSPACE/RadSun.csv";
        //String outPathToRadiationShadow			="output/GEOSPACE/RadShadow.csv";
        //String outPathToCanopy					="output/GEOSPACE/Canopy.csv";
        //String outPathToVPD						="output/GEOSPACE/VPD.csv";
		
        RichardsRootSolver1DMain R1DSolver     		= new RichardsRootSolver1DMain();
		GEOSPACEBuffer1D buffer 		  			= new GEOSPACEBuffer1D();
		WriteNetCDFGEOSPACE1DDouble writeNetCDF 	= new WriteNetCDFGEOSPACE1DDouble();
		ReadNetCDFGEOSPACEGrid1D readNetCDF    		= new ReadNetCDFGEOSPACEGrid1D();	
		
		InputReaderMain Input 											= new InputReaderMain();
		OutputWriterMain Output 										= new OutputWriterMain();
		JarvisNetRadiationStressFactorSolverMain JarvisStressFactor   	= new JarvisNetRadiationStressFactorSolverMain();
		PriestleyTaylorActualETSolverMain PtEt 							= new PriestleyTaylorActualETSolverMain();
		//PMEvaporationFromSoilAfterCanopySolverMain PMsoilevaporation 	= new PMEvaporationFromSoilAfterCanopySolverMain();
		//TotalEvapoTranspirationSolverMain TotalEvapoTranspiration 	= new TotalEvapoTranspirationSolverMain();
		InputDataMain InputBroker										= new InputDataMain();
		ETsBrokerOneFluxSolverMain ETsBrokerSolver 						= new ETsBrokerOneFluxSolverMain(); 
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
		
		
		Input.rootDepth  = -2;
		
		//JarvisStressFactor.etaE  = -0.5; //depth of the evaporation layer
		JarvisStressFactor.stressFactorModel = "LinearStressFactor";
		JarvisStressFactor.representativeStressFactorModel = "RootDensityWeightedMethod"; //SizeWightedMetod, AverageMetod
		//InputBroker.representativeEsModel = "";  	//SizeWaterWeightedMetod, AverageWaterWeightedMethod
		ETsBrokerSolver.representativeTsModel = "RootWaterWeightedMethod"; //SizeWaterWeightedMethod, AverageWaterWeightedMethod, RootWaterWeightedMethod //SizeWightedMethod, AverageMethod, RootWeightedMethod
       	
		
		Input.idCentroids="ID";
		Input.centroidElevation="Elevation";
		//Prospero.canopyHeight = 3.5;
		
		PtEt.alpha = 1.26;
        PtEt.soilFluxParameterDay = 0.35;
        PtEt.soilFluxParameterNight = 0.75;
		
		JarvisStressFactor.useRadiationStress    = false;
		JarvisStressFactor.useTemperatureStress  = false;
		JarvisStressFactor.useVDPStress		   = false;	
		JarvisStressFactor.useWaterStress = true;
		//Prospero.useEvaporationWaterStress = true;
		//JarvisStressFactor.defaultStress      = 1.0;
						
		JarvisStressFactor.alpha = 0.005;
		JarvisStressFactor.thetaR = 0.9;
		JarvisStressFactor.VPD0 = 5.0;
		JarvisStressFactor.Tl = -5.0;
		JarvisStressFactor.T0 = 20.0;
		JarvisStressFactor.Th = 45.0;
		//Prospero.typeOfCanopy="multilayer";

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
        //OmsTimeSeriesIteratorReader windReader 		 		= getTimeseriesReader(inPathToWind, fId, startDate, endDate, timeStepMinutes);
        //OmsTimeSeriesIteratorReader humidityReader 			= getTimeseriesReader(inPathToRelativeHumidity, fId, startDate, endDate, timeStepMinutes);
        //OmsTimeSeriesIteratorReader shortwaveReaderDirect 	= getTimeseriesReader(inPathToShortWaveRadiationDirect, fId, startDate, endDate,timeStepMinutes);
        //OmsTimeSeriesIteratorReader shortwaveReaderDiffuse 	= getTimeseriesReader(inPathToShortWaveRadiationDiffuse, fId, startDate, endDate,timeStepMinutes);
        //OmsTimeSeriesIteratorReader longwaveReader 			= getTimeseriesReader(inPathToLWRad, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader pressureReader 			= getTimeseriesReader(inPathToPressure, fId, startDate, endDate,timeStepMinutes);
        //OmsTimeSeriesIteratorReader leafAreaIndexReader		= getTimeseriesReader(inPathToLai, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader soilHeatFluxReader 		= getTimeseriesReader(inPathToSoilHeatFlux, fId, startDate, endDate,timeStepMinutes);
        OmsTimeSeriesIteratorReader netRadReader 			= getTimeseriesReader(inPathToNetRad, fId, startDate, endDate,timeStepMinutes);
		
		
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
		//InputBroker.etaE   			= JarvisStressFactor.etaE;
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
            Output.doPrintOutputPT = true;
            Input.tStartDate = startDate;
            Input.temporalStep = timeStepMinutes;

            
            
            soilHeatFluxReader.nextRecord();
            id2ValueMap = soilHeatFluxReader.outData;
            Input.inSoilFlux = id2ValueMap;
            
            pressureReader.nextRecord();
            id2ValueMap = pressureReader.outData;
            Input.inAtmosphericPressure = id2ValueMap;
            
            
            netRadReader.nextRecord();
            id2ValueMap = netRadReader.outData;
            Input.inNetRadiation = id2ValueMap;
            //Prospero.pm = pm;			

/////////////////////////////// SOLVE /////////////////////////
			
            R1DSolver.solve();

            JarvisStressFactor.theta = R1DSolver.thetasNew;
			
			Input.process();
			InputBroker.etaR = Input.defRootDepth;
			
			RootDensitySolver.solve();
			
			InputBroker.rootDensity = RootDensitySolver.defRootDensity;
		
			
			JarvisStressFactor.solve();
            
    
			//PMsoilevaporation.evaporationStressWater = JarvisStressFactor.evaporationStressWater;
			PtEt.stressFactor = JarvisStressFactor.stressSun;
            //Prospero.stressShade = JarvisStressFactor.stressShade;
            
			PtEt.process();
            
            //PMsoilevaporation.process();
            
            //TotalEvapoTranspiration.transpiration = Prospero.transpiration;
            //TotalEvapoTranspiration.evaporation = PMsoilevaporation.evaporation;
            //TotalEvapoTranspiration.process();
            
            InputBroker.g = JarvisStressFactor.g;
            InputBroker.GnT = JarvisStressFactor.GnT;
            //InputBroker.GnE = JarvisStressFactor.GnE;
            InputBroker.transpiration = PtEt.evapoTranspirationPT;
            //InputBroker.evaporation = PMsoilevaporation.evaporation;
           
            InputBroker.process();
			/////////////////////SOLO RICHARDS///////////////////////////////
			//ETsBrokerSolver.transpiration = 0;
			//ETsBrokerSolver.evaporation = 0;
            
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
			
			
            
            FluxEvapoTranspirationWriter.inData = Output.outLatentHeatPT;
            FluxEvapoTranspirationWriter.writeNextLine();	
			
			EvapoTranspirationWriter.inData = Output.outEvapoTranspirationPT;
			EvapoTranspirationWriter.writeNextLine();
			

			
	}	
	   		
		topBCReader.close();
		bottomBCReader.close();

		temperatureReader.close();        

        soilHeatFluxReader.close();
        pressureReader.close();

		FluxEvapoTranspirationWriter.close();


		EvapoTranspirationWriter.close();


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
