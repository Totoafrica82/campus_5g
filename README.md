# campus_5g
Campus 5G network proposal - qualitative assessment

To run the simulation, it is necessary to have the following MATLAB add-ons installed:
•	5G Toolbox
•	Antenna Toolbox
•	Communications Toolbox
•	Phased Array System Toolbox
•	Signal Processing Toolbox
•	Mapping Toolbox

3D OSM models from: https://osmbuildings.org/

The model **wat_coverage_final.m** incorporates multiple reflection orders to account for the complex multi-path environment created by campus buildings. Where environmental attenuation is integrated through a composite propagation model that accounts for atmospheric conditions relevant to the deployment area.

The simulation model **wat_compare_complex.m** was developed for the analysis and optimization of antenna configurations in our university's private 5G network deployment. The MATLAB based simulation platform enables quantitative evaluation of different antenna array architectures, beamforming techniques, and spatial multiplexing capabilities within the targeted coverage area. By using raytracing propagation models and detailed physical layer analysis, this simulation provides insights for optimal network deployment.

After each run the **wat_compare_complex.m** model saves antenna configuration in an additional file (eg. results_32TRX_128AE). After running all antenna configurations we move to the next simulation that compares all antenna systems **comparison_all.m**. 

This model visualizes the different antenna configurations advantages and helps with choosing the right configuration for our 5G network. The MIMO configuration analysis compares three antenna setups (8TRX, 32TRX, and 64TRX).

The **wat_KPI_complex_final.m **model simulates performance metrics for our 5G campus network implementation. This MATLAB simulation evaluates signal quality, throughput, latency, jitter, and reliability against established KPI targets across multiple antenna configurations. We strategically place UE (measuring points) within the campus environment to analyze performance and capacity under varying network loads. While the simulation generates detailed performance visualizations for all positioned UEs

Updated simulation **wat_KPI_1024qam_final.m** demonstrates the successful implementation of 1024QAM modulation into our private network model.


