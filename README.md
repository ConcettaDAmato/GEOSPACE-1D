# GEOSPACE-1D

GEOSPACE (Soil Plant Atmosphere Continuum Estimator model in GEOframe) is a Java-based ecohydrological modeling system developed within the GEOframe framework. It simulates the interactions between the atmosphere, vegetation, and soil, focusing on water and energy exchanges within the Earth’s Critical Zone. GEOSPACE integrates multiple models to provide a comprehensive understanding of these processes, including water flow, evapotranspiration (ET), and plant hydraulics.

This model adopts a **component-based** approach, breaking down the complexity of the **Soil-Plant-Atmosphere Continuum (SPAC)** into distinct, modular components. This approach enables the development of flexible, extensible, and reproducible models tailored to different user needs. Each component in GEOSPACE addresses a specific aspect of the SPAC, promoting interdisciplinary research and reliable simulations.

More details and information are available in [D’Amato, C.: Exploring the Soil-Plant-Atmosphere Continuum: Advancements, Integrated Modeling and Ecohydrological Insights, Ph.D. thesis, Center Agriculture Food Environment (C3A), University of Trento, 2024](https://abouthydrology.blogspot.com/2024/04/exploring-soil-plant-atmosphere.html).

## Scientific Overview

The SPAC system involves complex processes such as heat transfer, evapotranspiration, precipitation, water absorption, soil water flow, substance transport, and gas exchange. GEOSPACE captures these interactions using three main components:

1. **WHETGEO** (Water HEat Transport in GEOframe): Simulates water flow in the soil using the Richards-Richardson equation, solved with the Newton-Casulli-Zanolli algorithm, allowing for a robust treatment of variably saturated flow and heat transport [(Tubini and Rigon, 2022)](https://gmd.copernicus.org/articles/15/75/2022/gmd-15-75-2022.html).

2. **GEOET** (GEOframe EvapoTranspiration): A suite of evapotranspiration models including Priestley-Taylor, Penman-Monteith FAO, and Prospero, which incorporate energy budget calculations and detailed treatments of plant hydraulics and environmental stresses (e.g., Jarvis and Medlyn stomatal conductance models). 

3. **BrokerGEO**: A coupling component that manages data exchange between WHETGEO and GEOET, dynamically partitioning evapotranspiration, soil evaporation and plant transpiration across control volumes in the soil.

### Component-Based Design

GEOSPACE follows the principle of **modeling by components (MBC)** [(Rigon et al., 2022)](https://hess.copernicus.org/articles/26/4773/2022/). Each component is designed to perform a distinct task, allowing users to create customized modeling solutions by selecting and interchanging components as needed. This **object-oriented programming (OOP)** approach ensures that each component adheres to common interfaces, promoting flexibility and extensibility. The system also follows the "open to extensions, closed to modifications" principle, facilitating the addition of new features without altering existing components.

This modularity supports rigorous testing and inspection, ensuring that individual processes can be validated independently. Additionally, GEOSPACE supports the integration of new physical processes, such as alternative plant hydraulic models or new soil-root interaction formulations, making it a versatile tool for studying SPAC interactions under various environmental conditions.

### Key Features
- **Modular Architecture**: Each component (WHETGEO, GEOET, BrokerGEO) can be used independently or in combination, allowing users to adapt the system to their specific research needs.
- **Multiple ET Models**: GEOET offers a range of evapotranspiration models, from empirical to physically-based, enabling comparative studies of different formulations.
- **Dynamic Coupling**: BrokerGEO facilitates feedback between water flow and ET processes, ensuring that changes in soil moisture dynamically influence ET rates.
- **Environmental Stress Handling**: GEOET includes detailed models for environmental stresses (e.g., water scarcity) using functions like the Jarvis and Medlyn models.

### Applications
GEOSPACE is ideal for studying ecohydrological processes in the Earth’s Critical Zone. Applications include:
- **Water Resources Management**: Modeling water availability and plant-water relations under different environmental conditions.
- **Agriculture**: Simulating irrigation demands and the impact of water stress on crop yields.
- **Climate Change**: Investigating how changing climatic conditions affect the water cycle and plant growth.

## Software Structure
GEOSPACE is implemented in Java and integrated with the **OMS3** (Object Modeling System v3) framework, which supports modular environmental modeling. This framework allows users to customize and expand the model through the use of distinct components. GEOSPACE components are organized into Java packages that handle various hydrological, evapotranspiration, and data management tasks.

### GEOSPACE Components and Versions
GEOSPACE consists of multiple components, each versioned and packaged independently:

- **WHETGEO** (`whetgeo1d-1.2.9`): Simulates water and heat transport in the soil using the Richards equation for variably saturated flow.
  
- **BrokerGEO** (`brokergeo-1.3.9`): Manages the coupling between evapotranspiration and soil water content, partitioning ET into soil evaporation and transpiration.
  
- **GEOET** (`geoet-1.5.9`): Simulates evapotranspiration using various models (Priestley-Taylor, Penman-Monteith FAO, Prospero) and incorporates plant hydraulics and stress factor calculations.
  
- **Buffer** (`buffer-1.1.9`): Manages data in memory and handles printing tasks during simulations.
  
- **ClosureEquation** (`closureequation-1.1.9`): Manages hydraulic and thermal properties equations used by WHETGEO.
  
- **NetCDF** (`netcdf-1.1.9`): Handles input and output operations for NetCDF file formats, facilitating large dataset management.
  
- **Numerical** (`numerical-1.0.2`): Provides core numerical algorithms for solving linear and nonlinear systems.

The system supports a high degree of modularity, allowing new processes or models to be integrated easily by creating new components or extending existing ones. The architecture also ensures that simulations are reproducible, as OMS3 workflows are stored in `.sim` files that preserve all the configuration details for exact replication.

## Source Code
- [WHETGEO Repository](https://github.com/geoframecomponents/WHETGEO-1D)
- [GEOET Repository](https://github.com/geoframecomponents/GEOET)
- [BrokerGEO Repository](https://github.com/geoframecomponents/BrokerGEO)

## Executables
- [GEOSPACE Executables](https://github.com/GEOframeOMSProjects/OMS_Project_GEOSPACE-1D)
- [WHETGEO Executables](https://github.com/GEOframeOMSProjects/OMS_Project_WHETGEO1D)
- [GEOET Executables](https://github.com/GEOframeOMSProjects/OMS_Project_GEOET)

## Contributing
Contributions are welcome! Please adhere to the [GEOframe community guidelines](https://geoframe.blogspot.com/2020/05/geoframe-community-publication-policy.html) when submitting your contributions.

## Acknowledgements
- **Concetta D’Amato**, **Niccolò Tubini**, and **Riccardo Rigon** contributed significantly to the development of GEOSPACE and its components.
- GEOSPACE development was supported by Concetta D'Amato Ph.D. grant by [C3A-UniTrento](https://www.centro3a.unitn.it/) and by WATZON, WATERSTEM and WATSON projects. 
  
