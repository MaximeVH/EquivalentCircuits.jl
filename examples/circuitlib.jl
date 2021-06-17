#Running this code will take a while, it is included here for illustration purposed. The resulting circuit library CSV file (Circuitlibrary.csv) can already conveniently
#be found in the root folder.

initialize_circuitlibrary("R1-[C2,R3-[C4,R5]]",[20,4e-9,3400,4e-6,2500],"Nonbiological","https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/")
add_to_circuitlibrary("[R1,[C2,R3]-R4]",[4310,3.5e-9,2.9e6,220],"Biological","Electrical Impedance Studies on Potato and Alfalfa Tissue, Hayden et al")
add_to_circuitlibrary("[R1,C2-[R3,C4-R5]]",[7355,5.11e-9,3274,5e-9,346],"Biological","Electrical Impedance Analysis in Plant Tissues: A Double Shell Model, Zhang et al")
add_to_circuitlibrary("R1-[C2,R3]-[C4,R5]","Biological","Bioimpedance analysis of vascular tissue and fluid flow in human and plant body: A review, Prasad et al (plant cross-section)")
add_to_circuitlibrary("[R1-C2,R3]","Biological","Bioimpedance analysis of vascular tissue and fluid flow in human and plant body: A review, Prasad et al (body segment tissue model)")
add_to_circuitlibrary("R1-L2-[C3,R4]","Biological","Bioimpedance analysis of vascular tissue and fluid flow in human and plant body: A review, Prasad et al (arterial system RLCR model)")
add_to_circuitlibrary("P1-R2-[C3,R4]","Biological","Electrochemical impedance spectroscopy as a versatile tool for the characterization of neural tissue: A mini review, Krukiewicz et al")
add_to_circuitlibrary("[R1,C2-R3-[C4,R5]]",[2.6e3,1e-6,470,10,10e3],"Biological","Impedance spectroscopy on living tissue for determination of the state of organs, Gersing et al (porcine liver)")
add_to_circuitlibrary("[R1-[P2,R3],C4]","Biological","A suggested circuit for bacterial biofilms (Erkuden)")
add_to_circuitlibrary("R1-[R2,C3]",[20,250,40e-6],"Nonbiological","Randles circuit, https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/")
add_to_circuitlibrary("R1-[C2,R3-P4]",[20,40e-6,250,[150,0.5]],"Nonbiological"," Mixed Kinetic and Charge-Transfer Control Randles cell, https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/")
add_to_circuitlibrary("R1-C2-[R3,P4]","Biological"," Rapid detection of bacterial proliferation in food samples,Puttaswamy et al")
add_to_circuitlibrary("R1-[P2,R3]",[234,[2e8,0.56],15789],"Biological"," Electrical impedance spectroscopic study of mandarin orange during ripening, Chowdhury et al")
add_encoding_to_circuitlibrary("+++--+PPRPRLR","Nonbiological","Characterization of high-power lithium-ion batteries by electrochemical impedance spectroscopy, Andre et al")
add_to_circuitlibrary("[C1,R2-C3]",[2.1e-12,1.6e3,2.2e-12],"Nonbiological","Impedance spectroscopy study of hardened Portland cement paste, cabeza et al")
