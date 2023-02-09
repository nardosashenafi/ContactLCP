  function targMap = targDataMap(),

  ;%***********************
  ;% Create Parameter Map *
  ;%***********************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 7;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc paramMap
    ;%
    paramMap.nSections           = nTotSects;
    paramMap.sectIdxOffset       = sectIdxOffset;
      paramMap.sections(nTotSects) = dumSection; %prealloc
    paramMap.nTotData            = -1;
    
    ;%
    ;% Auto data (parameter_estimation_hardware_P)
    ;%
      section.nData     = 17;
      section.data(17)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.A0
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.Kt
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardware_P.Nenc
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.R
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardware_P.T
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardware_P.Vmax
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardware_P.af
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% parameter_estimation_hardware_P.bf
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 9;
	
	  ;% parameter_estimation_hardware_P.bhat0
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 12;
	
	  ;% parameter_estimation_hardware_P.k
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 13;
	
	  ;% parameter_estimation_hardware_P.lambda
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 14;
	
	  ;% parameter_estimation_hardware_P.mhat0
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 15;
	
	  ;% parameter_estimation_hardware_P.omega
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 16;
	
	  ;% parameter_estimation_hardware_P.phi
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 17;
	
	  ;% parameter_estimation_hardware_P.r_enc
	  section.data(15).logicalSrcIdx = 14;
	  section.data(15).dtTransOffset = 18;
	
	  ;% parameter_estimation_hardware_P.rg
	  section.data(16).logicalSrcIdx = 15;
	  section.data(16).dtTransOffset = 19;
	
	  ;% parameter_estimation_hardware_P.rm
	  section.data(17).logicalSrcIdx = 16;
	  section.data(17).dtTransOffset = 20;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILReadEncoderTimebase1_clock
	  section.data(1).logicalSrcIdx = 17;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 3;
      section.data(3)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILReadEncoderTimebase1_channel
	  section.data(1).logicalSrcIdx = 18;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.HILWriteAnalog_channels
	  section.data(2).logicalSrcIdx = 19;
	  section.data(2).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.HILReadEncoderTimebase1_samples
	  section.data(3).logicalSrcIdx = 20;
	  section.data(3).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 25;
      section.data(25)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILInitialize_OOTerminate
	  section.data(1).logicalSrcIdx = 21;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_OOExit
	  section.data(2).logicalSrcIdx = 22;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_OOStart
	  section.data(3).logicalSrcIdx = 23;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_OOEnter
	  section.data(4).logicalSrcIdx = 24;
	  section.data(4).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOFinal
	  section.data(5).logicalSrcIdx = 25;
	  section.data(5).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POFinal
	  section.data(6).logicalSrcIdx = 26;
	  section.data(6).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AIHigh
	  section.data(7).logicalSrcIdx = 27;
	  section.data(7).dtTransOffset = 6;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AILow
	  section.data(8).logicalSrcIdx = 28;
	  section.data(8).dtTransOffset = 7;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOHigh
	  section.data(9).logicalSrcIdx = 29;
	  section.data(9).dtTransOffset = 8;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOLow
	  section.data(10).logicalSrcIdx = 30;
	  section.data(10).dtTransOffset = 9;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOInitial
	  section.data(11).logicalSrcIdx = 31;
	  section.data(11).dtTransOffset = 10;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOWatchdog
	  section.data(12).logicalSrcIdx = 32;
	  section.data(12).dtTransOffset = 11;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POFrequency
	  section.data(13).logicalSrcIdx = 33;
	  section.data(13).dtTransOffset = 12;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POInitial
	  section.data(14).logicalSrcIdx = 34;
	  section.data(14).dtTransOffset = 13;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POWatchdog
	  section.data(15).logicalSrcIdx = 35;
	  section.data(15).dtTransOffset = 14;
	
	  ;% parameter_estimation_hardware_P.DiscreteFilter_NumCoef
	  section.data(16).logicalSrcIdx = 36;
	  section.data(16).dtTransOffset = 15;
	
	  ;% parameter_estimation_hardware_P.DiscreteFilter_DenCoef
	  section.data(17).logicalSrcIdx = 37;
	  section.data(17).dtTransOffset = 17;
	
	  ;% parameter_estimation_hardware_P.DiscreteFilter_InitialStates
	  section.data(18).logicalSrcIdx = 38;
	  section.data(18).dtTransOffset = 19;
	
	  ;% parameter_estimation_hardware_P.DiscreteFilter1_InitialStates
	  section.data(19).logicalSrcIdx = 39;
	  section.data(19).dtTransOffset = 20;
	
	  ;% parameter_estimation_hardware_P.SineWave2_Bias
	  section.data(20).logicalSrcIdx = 40;
	  section.data(20).dtTransOffset = 21;
	
	  ;% parameter_estimation_hardware_P.SineWave1_Bias
	  section.data(21).logicalSrcIdx = 41;
	  section.data(21).dtTransOffset = 22;
	
	  ;% parameter_estimation_hardware_P.SineWave_Bias
	  section.data(22).logicalSrcIdx = 42;
	  section.data(22).dtTransOffset = 23;
	
	  ;% parameter_estimation_hardware_P.Onoroff_Gain
	  section.data(23).logicalSrcIdx = 43;
	  section.data(23).dtTransOffset = 24;
	
	  ;% parameter_estimation_hardware_P.Gain1_Gain
	  section.data(24).logicalSrcIdx = 44;
	  section.data(24).dtTransOffset = 25;
	
	  ;% parameter_estimation_hardware_P.Gain3_Gain
	  section.data(25).logicalSrcIdx = 45;
	  section.data(25).dtTransOffset = 26;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILInitialize_CKChannels
	  section.data(1).logicalSrcIdx = 46;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOWatchdog
	  section.data(2).logicalSrcIdx = 47;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIInitial
	  section.data(3).logicalSrcIdx = 48;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POModes
	  section.data(4).logicalSrcIdx = 49;
	  section.data(4).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILInitialize_AIChannels
	  section.data(1).logicalSrcIdx = 50;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOChannels
	  section.data(2).logicalSrcIdx = 51;
	  section.data(2).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIChannels
	  section.data(3).logicalSrcIdx = 52;
	  section.data(3).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIQuadrature
	  section.data(4).logicalSrcIdx = 53;
	  section.data(4).dtTransOffset = 6;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
      clear section
      
      section.nData     = 37;
      section.data(37)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_P.HILInitialize_Active
	  section.data(1).logicalSrcIdx = 54;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOTerminate
	  section.data(2).logicalSrcIdx = 55;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOExit
	  section.data(3).logicalSrcIdx = 56;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOTerminate
	  section.data(4).logicalSrcIdx = 57;
	  section.data(4).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOExit
	  section.data(5).logicalSrcIdx = 58;
	  section.data(5).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POTerminate
	  section.data(6).logicalSrcIdx = 59;
	  section.data(6).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POExit
	  section.data(7).logicalSrcIdx = 60;
	  section.data(7).dtTransOffset = 6;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_CKPStart
	  section.data(8).logicalSrcIdx = 61;
	  section.data(8).dtTransOffset = 7;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_CKPEnter
	  section.data(9).logicalSrcIdx = 62;
	  section.data(9).dtTransOffset = 8;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_CKStart
	  section.data(10).logicalSrcIdx = 63;
	  section.data(10).dtTransOffset = 9;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_CKEnter
	  section.data(11).logicalSrcIdx = 64;
	  section.data(11).dtTransOffset = 10;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AIPStart
	  section.data(12).logicalSrcIdx = 65;
	  section.data(12).dtTransOffset = 11;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AIPEnter
	  section.data(13).logicalSrcIdx = 66;
	  section.data(13).dtTransOffset = 12;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOPStart
	  section.data(14).logicalSrcIdx = 67;
	  section.data(14).dtTransOffset = 13;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOPEnter
	  section.data(15).logicalSrcIdx = 68;
	  section.data(15).dtTransOffset = 14;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOStart
	  section.data(16).logicalSrcIdx = 69;
	  section.data(16).dtTransOffset = 15;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOEnter
	  section.data(17).logicalSrcIdx = 70;
	  section.data(17).dtTransOffset = 16;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_AOReset
	  section.data(18).logicalSrcIdx = 71;
	  section.data(18).dtTransOffset = 17;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOPStart
	  section.data(19).logicalSrcIdx = 72;
	  section.data(19).dtTransOffset = 18;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOPEnter
	  section.data(20).logicalSrcIdx = 73;
	  section.data(20).dtTransOffset = 19;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOStart
	  section.data(21).logicalSrcIdx = 74;
	  section.data(21).dtTransOffset = 20;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOEnter
	  section.data(22).logicalSrcIdx = 75;
	  section.data(22).dtTransOffset = 21;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOReset
	  section.data(23).logicalSrcIdx = 76;
	  section.data(23).dtTransOffset = 22;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIPStart
	  section.data(24).logicalSrcIdx = 77;
	  section.data(24).dtTransOffset = 23;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIPEnter
	  section.data(25).logicalSrcIdx = 78;
	  section.data(25).dtTransOffset = 24;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIStart
	  section.data(26).logicalSrcIdx = 79;
	  section.data(26).dtTransOffset = 25;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_EIEnter
	  section.data(27).logicalSrcIdx = 80;
	  section.data(27).dtTransOffset = 26;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POPStart
	  section.data(28).logicalSrcIdx = 81;
	  section.data(28).dtTransOffset = 27;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POPEnter
	  section.data(29).logicalSrcIdx = 82;
	  section.data(29).dtTransOffset = 28;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POStart
	  section.data(30).logicalSrcIdx = 83;
	  section.data(30).dtTransOffset = 29;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POEnter
	  section.data(31).logicalSrcIdx = 84;
	  section.data(31).dtTransOffset = 30;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_POReset
	  section.data(32).logicalSrcIdx = 85;
	  section.data(32).dtTransOffset = 31;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_OOReset
	  section.data(33).logicalSrcIdx = 86;
	  section.data(33).dtTransOffset = 32;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOFinal
	  section.data(34).logicalSrcIdx = 87;
	  section.data(34).dtTransOffset = 33;
	
	  ;% parameter_estimation_hardware_P.HILInitialize_DOInitial
	  section.data(35).logicalSrcIdx = 88;
	  section.data(35).dtTransOffset = 34;
	
	  ;% parameter_estimation_hardware_P.HILReadEncoderTimebase1_Active
	  section.data(36).logicalSrcIdx = 89;
	  section.data(36).dtTransOffset = 35;
	
	  ;% parameter_estimation_hardware_P.HILWriteAnalog_Active
	  section.data(37).logicalSrcIdx = 90;
	  section.data(37).dtTransOffset = 36;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(7) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (parameter)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    paramMap.nTotData = nTotData;
    


  ;%**************************
  ;% Create Block Output Map *
  ;%**************************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 1;
    sectIdxOffset = 0;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc sigMap
    ;%
    sigMap.nSections           = nTotSects;
    sigMap.sectIdxOffset       = sectIdxOffset;
      sigMap.sections(nTotSects) = dumSection; %prealloc
    sigMap.nTotData            = -1;
    
    ;%
    ;% Auto data (parameter_estimation_hardware_B)
    ;%
      section.nData     = 11;
      section.data(11)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardware_B.Kenc
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardware_B.DiscreteFilter1
	  section.data(2).logicalSrcIdx = 2;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardware_B.Integrator
	  section.data(3).logicalSrcIdx = 3;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardware_B.Add1
	  section.data(4).logicalSrcIdx = 4;
	  section.data(4).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardware_B.Integrator1
	  section.data(5).logicalSrcIdx = 5;
	  section.data(5).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardware_B.Add
	  section.data(6).logicalSrcIdx = 6;
	  section.data(6).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardware_B.Gain6
	  section.data(7).logicalSrcIdx = 7;
	  section.data(7).dtTransOffset = 6;
	
	  ;% parameter_estimation_hardware_B.Saturation1
	  section.data(8).logicalSrcIdx = 8;
	  section.data(8).dtTransOffset = 7;
	
	  ;% parameter_estimation_hardware_B.Onoroff
	  section.data(9).logicalSrcIdx = 9;
	  section.data(9).dtTransOffset = 8;
	
	  ;% parameter_estimation_hardware_B.Product
	  section.data(10).logicalSrcIdx = 10;
	  section.data(10).dtTransOffset = 9;
	
	  ;% parameter_estimation_hardware_B.Product1
	  section.data(11).logicalSrcIdx = 11;
	  section.data(11).dtTransOffset = 10;
	
      nTotData = nTotData + section.nData;
      sigMap.sections(1) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (signal)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    sigMap.nTotData = nTotData;
    


  ;%*******************
  ;% Create DWork Map *
  ;%*******************
      
    nTotData      = 0; %add to this count as we go
    nTotSects     = 5;
    sectIdxOffset = 1;
    
    ;%
    ;% Define dummy sections & preallocate arrays
    ;%
    dumSection.nData = -1;  
    dumSection.data  = [];
    
    dumData.logicalSrcIdx = -1;
    dumData.dtTransOffset = -1;
    
    ;%
    ;% Init/prealloc dworkMap
    ;%
    dworkMap.nSections           = nTotSects;
    dworkMap.sectIdxOffset       = sectIdxOffset;
      dworkMap.sections(nTotSects) = dumSection; %prealloc
    dworkMap.nTotData            = -1;
    
    ;%
    ;% Auto data (parameter_estimation_hardwar_DW)
    ;%
      section.nData     = 10;
      section.data(10)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardwar_DW.DiscreteFilter_states
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardwar_DW.DiscreteFilter1_states
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_AIMinimums
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_AIMaximums
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_AOMinimums
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 7;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_AOMaximums
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 9;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_AOVoltages
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 11;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_FilterFrequency
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 13;
	
	  ;% parameter_estimation_hardwar_DW.DiscreteFilter_tmp
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 15;
	
	  ;% parameter_estimation_hardwar_DW.DiscreteFilter1_tmp
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 16;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardwar_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 10;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardwar_DW.HILReadEncoderTimebase1_Task
	  section.data(1).logicalSrcIdx = 11;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 9;
      section.data(9)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardwar_DW.HILWriteAnalog_PWORK
	  section.data(1).logicalSrcIdx = 12;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardwar_DW.Voltage_PWORK.LoggedData
	  section.data(2).logicalSrcIdx = 13;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardwar_DW.bhat_PWORK.LoggedData
	  section.data(3).logicalSrcIdx = 14;
	  section.data(3).dtTransOffset = 2;
	
	  ;% parameter_estimation_hardwar_DW.mhat_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 15;
	  section.data(4).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardwar_DW.v_PWORK.LoggedData
	  section.data(5).logicalSrcIdx = 16;
	  section.data(5).dtTransOffset = 4;
	
	  ;% parameter_estimation_hardwar_DW.voltage_PWORK.LoggedData
	  section.data(6).logicalSrcIdx = 17;
	  section.data(6).dtTransOffset = 5;
	
	  ;% parameter_estimation_hardwar_DW.x_PWORK.LoggedData
	  section.data(7).logicalSrcIdx = 18;
	  section.data(7).dtTransOffset = 6;
	
	  ;% parameter_estimation_hardwar_DW.xtilde_PWORK.LoggedData
	  section.data(8).logicalSrcIdx = 19;
	  section.data(8).dtTransOffset = 7;
	
	  ;% parameter_estimation_hardwar_DW.xtildedot_PWORK.LoggedData
	  section.data(9).logicalSrcIdx = 20;
	  section.data(9).dtTransOffset = 8;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% parameter_estimation_hardwar_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 21;
	  section.data(1).dtTransOffset = 0;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_QuadratureModes
	  section.data(2).logicalSrcIdx = 22;
	  section.data(2).dtTransOffset = 1;
	
	  ;% parameter_estimation_hardwar_DW.HILInitialize_InitialEICounts
	  section.data(3).logicalSrcIdx = 23;
	  section.data(3).dtTransOffset = 3;
	
	  ;% parameter_estimation_hardwar_DW.HILReadEncoderTimebase1_Buffer
	  section.data(4).logicalSrcIdx = 24;
	  section.data(4).dtTransOffset = 5;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
    
      ;%
      ;% Non-auto Data (dwork)
      ;%
    

    ;%
    ;% Add final counts to struct.
    ;%
    dworkMap.nTotData = nTotData;
    


  ;%
  ;% Add individual maps to base struct.
  ;%

  targMap.paramMap  = paramMap;    
  targMap.signalMap = sigMap;
  targMap.dworkMap  = dworkMap;
  
  ;%
  ;% Add checksums to base struct.
  ;%


  targMap.checksum0 = 2267798225;
  targMap.checksum1 = 1350053079;
  targMap.checksum2 = 1186518877;
  targMap.checksum3 = 3431875335;

