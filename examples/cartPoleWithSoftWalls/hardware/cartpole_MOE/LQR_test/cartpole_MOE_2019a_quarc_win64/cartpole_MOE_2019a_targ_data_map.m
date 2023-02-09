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
    ;% Auto data (cartpole_MOE_2019a_P)
    ;%
      section.nData     = 11;
      section.data(11)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.Fmax
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.Kt
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_P.Nenc
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.R
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_P.T
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_P.Vmax
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 5;
	
	  ;% cartpole_MOE_2019a_P.af
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_P.bf
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 9;
	
	  ;% cartpole_MOE_2019a_P.r_enc
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 12;
	
	  ;% cartpole_MOE_2019a_P.rg
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 13;
	
	  ;% cartpole_MOE_2019a_P.rm
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 14;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILReadEncoderTimebase_clock
	  section.data(1).logicalSrcIdx = 11;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(2) = section;
      clear section
      
      section.nData     = 3;
      section.data(3)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILReadEncoderTimebase_channels
	  section.data(1).logicalSrcIdx = 12;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.HILWriteAnalog_channels
	  section.data(2).logicalSrcIdx = 13;
	  section.data(2).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.HILReadEncoderTimebase_samples_
	  section.data(3).logicalSrcIdx = 14;
	  section.data(3).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(3) = section;
      clear section
      
      section.nData     = 38;
      section.data(38)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILInitialize_OOTerminate
	  section.data(1).logicalSrcIdx = 15;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_OOExit
	  section.data(2).logicalSrcIdx = 16;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_OOStart
	  section.data(3).logicalSrcIdx = 17;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_OOEnter
	  section.data(4).logicalSrcIdx = 18;
	  section.data(4).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOFinal
	  section.data(5).logicalSrcIdx = 19;
	  section.data(5).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POFinal
	  section.data(6).logicalSrcIdx = 20;
	  section.data(6).dtTransOffset = 5;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AIHigh
	  section.data(7).logicalSrcIdx = 21;
	  section.data(7).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AILow
	  section.data(8).logicalSrcIdx = 22;
	  section.data(8).dtTransOffset = 7;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOHigh
	  section.data(9).logicalSrcIdx = 23;
	  section.data(9).dtTransOffset = 8;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOLow
	  section.data(10).logicalSrcIdx = 24;
	  section.data(10).dtTransOffset = 9;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOInitial
	  section.data(11).logicalSrcIdx = 25;
	  section.data(11).dtTransOffset = 10;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOWatchdog
	  section.data(12).logicalSrcIdx = 26;
	  section.data(12).dtTransOffset = 11;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POFrequency
	  section.data(13).logicalSrcIdx = 27;
	  section.data(13).dtTransOffset = 12;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POInitial
	  section.data(14).logicalSrcIdx = 28;
	  section.data(14).dtTransOffset = 13;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POWatchdog
	  section.data(15).logicalSrcIdx = 29;
	  section.data(15).dtTransOffset = 14;
	
	  ;% cartpole_MOE_2019a_P.Gain4_Gain
	  section.data(16).logicalSrcIdx = 30;
	  section.data(16).dtTransOffset = 15;
	
	  ;% cartpole_MOE_2019a_P.Constant_Value
	  section.data(17).logicalSrcIdx = 31;
	  section.data(17).dtTransOffset = 16;
	
	  ;% cartpole_MOE_2019a_P.Gain5_Gain
	  section.data(18).logicalSrcIdx = 32;
	  section.data(18).dtTransOffset = 17;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_NumCoef
	  section.data(19).logicalSrcIdx = 33;
	  section.data(19).dtTransOffset = 18;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_DenCoef
	  section.data(20).logicalSrcIdx = 34;
	  section.data(20).dtTransOffset = 20;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_InitialStates
	  section.data(21).logicalSrcIdx = 35;
	  section.data(21).dtTransOffset = 22;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter1_InitialStates
	  section.data(22).logicalSrcIdx = 36;
	  section.data(22).dtTransOffset = 23;
	
	  ;% cartpole_MOE_2019a_P.Gain6_Gain
	  section.data(23).logicalSrcIdx = 37;
	  section.data(23).dtTransOffset = 24;
	
	  ;% cartpole_MOE_2019a_P.Constant1_Value
	  section.data(24).logicalSrcIdx = 38;
	  section.data(24).dtTransOffset = 25;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_NumCoef_b
	  section.data(25).logicalSrcIdx = 39;
	  section.data(25).dtTransOffset = 26;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_DenCoef_f
	  section.data(26).logicalSrcIdx = 40;
	  section.data(26).dtTransOffset = 28;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_InitialStates_g
	  section.data(27).logicalSrcIdx = 41;
	  section.data(27).dtTransOffset = 30;
	
	  ;% cartpole_MOE_2019a_P.DiscreteFilter_InitialStates_h
	  section.data(28).logicalSrcIdx = 42;
	  section.data(28).dtTransOffset = 31;
	
	  ;% cartpole_MOE_2019a_P.Gain7_Gain
	  section.data(29).logicalSrcIdx = 43;
	  section.data(29).dtTransOffset = 32;
	
	  ;% cartpole_MOE_2019a_P.Gain_Gain
	  section.data(30).logicalSrcIdx = 44;
	  section.data(30).dtTransOffset = 33;
	
	  ;% cartpole_MOE_2019a_P.Integrator_IC
	  section.data(31).logicalSrcIdx = 45;
	  section.data(31).dtTransOffset = 34;
	
	  ;% cartpole_MOE_2019a_P.Gain1_Gain
	  section.data(32).logicalSrcIdx = 46;
	  section.data(32).dtTransOffset = 35;
	
	  ;% cartpole_MOE_2019a_P.Saturation2_UpperSat
	  section.data(33).logicalSrcIdx = 47;
	  section.data(33).dtTransOffset = 36;
	
	  ;% cartpole_MOE_2019a_P.Saturation2_LowerSat
	  section.data(34).logicalSrcIdx = 48;
	  section.data(34).dtTransOffset = 37;
	
	  ;% cartpole_MOE_2019a_P.Onoroff1_Gain
	  section.data(35).logicalSrcIdx = 49;
	  section.data(35).dtTransOffset = 38;
	
	  ;% cartpole_MOE_2019a_P.Onoroff_Gain
	  section.data(36).logicalSrcIdx = 50;
	  section.data(36).dtTransOffset = 39;
	
	  ;% cartpole_MOE_2019a_P.Gain3_Gain
	  section.data(37).logicalSrcIdx = 51;
	  section.data(37).dtTransOffset = 40;
	
	  ;% cartpole_MOE_2019a_P.Gain2_Gain
	  section.data(38).logicalSrcIdx = 52;
	  section.data(38).dtTransOffset = 41;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(4) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILInitialize_CKChannels
	  section.data(1).logicalSrcIdx = 53;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOWatchdog
	  section.data(2).logicalSrcIdx = 54;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIInitial
	  section.data(3).logicalSrcIdx = 55;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POModes
	  section.data(4).logicalSrcIdx = 56;
	  section.data(4).dtTransOffset = 3;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(5) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILInitialize_AIChannels
	  section.data(1).logicalSrcIdx = 57;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOChannels
	  section.data(2).logicalSrcIdx = 58;
	  section.data(2).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIChannels
	  section.data(3).logicalSrcIdx = 59;
	  section.data(3).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIQuadrature
	  section.data(4).logicalSrcIdx = 60;
	  section.data(4).dtTransOffset = 6;
	
      nTotData = nTotData + section.nData;
      paramMap.sections(6) = section;
      clear section
      
      section.nData     = 37;
      section.data(37)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_P.HILInitialize_Active
	  section.data(1).logicalSrcIdx = 61;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOTerminate
	  section.data(2).logicalSrcIdx = 62;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOExit
	  section.data(3).logicalSrcIdx = 63;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOTerminate
	  section.data(4).logicalSrcIdx = 64;
	  section.data(4).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOExit
	  section.data(5).logicalSrcIdx = 65;
	  section.data(5).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POTerminate
	  section.data(6).logicalSrcIdx = 66;
	  section.data(6).dtTransOffset = 5;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POExit
	  section.data(7).logicalSrcIdx = 67;
	  section.data(7).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_CKPStart
	  section.data(8).logicalSrcIdx = 68;
	  section.data(8).dtTransOffset = 7;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_CKPEnter
	  section.data(9).logicalSrcIdx = 69;
	  section.data(9).dtTransOffset = 8;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_CKStart
	  section.data(10).logicalSrcIdx = 70;
	  section.data(10).dtTransOffset = 9;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_CKEnter
	  section.data(11).logicalSrcIdx = 71;
	  section.data(11).dtTransOffset = 10;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AIPStart
	  section.data(12).logicalSrcIdx = 72;
	  section.data(12).dtTransOffset = 11;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AIPEnter
	  section.data(13).logicalSrcIdx = 73;
	  section.data(13).dtTransOffset = 12;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOPStart
	  section.data(14).logicalSrcIdx = 74;
	  section.data(14).dtTransOffset = 13;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOPEnter
	  section.data(15).logicalSrcIdx = 75;
	  section.data(15).dtTransOffset = 14;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOStart
	  section.data(16).logicalSrcIdx = 76;
	  section.data(16).dtTransOffset = 15;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOEnter
	  section.data(17).logicalSrcIdx = 77;
	  section.data(17).dtTransOffset = 16;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_AOReset
	  section.data(18).logicalSrcIdx = 78;
	  section.data(18).dtTransOffset = 17;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOPStart
	  section.data(19).logicalSrcIdx = 79;
	  section.data(19).dtTransOffset = 18;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOPEnter
	  section.data(20).logicalSrcIdx = 80;
	  section.data(20).dtTransOffset = 19;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOStart
	  section.data(21).logicalSrcIdx = 81;
	  section.data(21).dtTransOffset = 20;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOEnter
	  section.data(22).logicalSrcIdx = 82;
	  section.data(22).dtTransOffset = 21;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOReset
	  section.data(23).logicalSrcIdx = 83;
	  section.data(23).dtTransOffset = 22;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIPStart
	  section.data(24).logicalSrcIdx = 84;
	  section.data(24).dtTransOffset = 23;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIPEnter
	  section.data(25).logicalSrcIdx = 85;
	  section.data(25).dtTransOffset = 24;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIStart
	  section.data(26).logicalSrcIdx = 86;
	  section.data(26).dtTransOffset = 25;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_EIEnter
	  section.data(27).logicalSrcIdx = 87;
	  section.data(27).dtTransOffset = 26;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POPStart
	  section.data(28).logicalSrcIdx = 88;
	  section.data(28).dtTransOffset = 27;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POPEnter
	  section.data(29).logicalSrcIdx = 89;
	  section.data(29).dtTransOffset = 28;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POStart
	  section.data(30).logicalSrcIdx = 90;
	  section.data(30).dtTransOffset = 29;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POEnter
	  section.data(31).logicalSrcIdx = 91;
	  section.data(31).dtTransOffset = 30;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_POReset
	  section.data(32).logicalSrcIdx = 92;
	  section.data(32).dtTransOffset = 31;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_OOReset
	  section.data(33).logicalSrcIdx = 93;
	  section.data(33).dtTransOffset = 32;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOFinal
	  section.data(34).logicalSrcIdx = 94;
	  section.data(34).dtTransOffset = 33;
	
	  ;% cartpole_MOE_2019a_P.HILInitialize_DOInitial
	  section.data(35).logicalSrcIdx = 95;
	  section.data(35).dtTransOffset = 34;
	
	  ;% cartpole_MOE_2019a_P.HILReadEncoderTimebase_Active
	  section.data(36).logicalSrcIdx = 96;
	  section.data(36).dtTransOffset = 35;
	
	  ;% cartpole_MOE_2019a_P.HILWriteAnalog_Active
	  section.data(37).logicalSrcIdx = 97;
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
    ;% Auto data (cartpole_MOE_2019a_B)
    ;%
      section.nData     = 14;
      section.data(14)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_B.Kenc
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_B.DiscreteFilter1
	  section.data(2).logicalSrcIdx = 2;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_B.Sum1
	  section.data(3).logicalSrcIdx = 3;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_B.DiscreteFilter
	  section.data(4).logicalSrcIdx = 5;
	  section.data(4).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_B.Gain
	  section.data(5).logicalSrcIdx = 6;
	  section.data(5).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_B.Saturation2
	  section.data(6).logicalSrcIdx = 7;
	  section.data(6).dtTransOffset = 5;
	
	  ;% cartpole_MOE_2019a_B.Saturation1
	  section.data(7).logicalSrcIdx = 8;
	  section.data(7).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_B.Onoroff1
	  section.data(8).logicalSrcIdx = 9;
	  section.data(8).dtTransOffset = 7;
	
	  ;% cartpole_MOE_2019a_B.Gain6
	  section.data(9).logicalSrcIdx = 10;
	  section.data(9).dtTransOffset = 8;
	
	  ;% cartpole_MOE_2019a_B.Saturation1_k
	  section.data(10).logicalSrcIdx = 11;
	  section.data(10).dtTransOffset = 9;
	
	  ;% cartpole_MOE_2019a_B.Onoroff
	  section.data(11).logicalSrcIdx = 12;
	  section.data(11).dtTransOffset = 10;
	
	  ;% cartpole_MOE_2019a_B.u60Nenc
	  section.data(12).logicalSrcIdx = 13;
	  section.data(12).dtTransOffset = 11;
	
	  ;% cartpole_MOE_2019a_B.Sum1_o
	  section.data(13).logicalSrcIdx = 14;
	  section.data(13).dtTransOffset = 12;
	
	  ;% cartpole_MOE_2019a_B.Sum
	  section.data(14).logicalSrcIdx = 15;
	  section.data(14).dtTransOffset = 13;
	
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
    nTotSects     = 6;
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
    ;% Auto data (cartpole_MOE_2019a_DW)
    ;%
      section.nData     = 14;
      section.data(14)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_states
	  section.data(1).logicalSrcIdx = 0;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter1_states
	  section.data(2).logicalSrcIdx = 1;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_states_p
	  section.data(3).logicalSrcIdx = 2;
	  section.data(3).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_states_a
	  section.data(4).logicalSrcIdx = 3;
	  section.data(4).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_AIMinimums
	  section.data(5).logicalSrcIdx = 4;
	  section.data(5).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_AIMaximums
	  section.data(6).logicalSrcIdx = 5;
	  section.data(6).dtTransOffset = 8;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_AOMinimums
	  section.data(7).logicalSrcIdx = 6;
	  section.data(7).dtTransOffset = 10;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_AOMaximums
	  section.data(8).logicalSrcIdx = 7;
	  section.data(8).dtTransOffset = 12;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_AOVoltages
	  section.data(9).logicalSrcIdx = 8;
	  section.data(9).dtTransOffset = 14;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_FilterFrequency
	  section.data(10).logicalSrcIdx = 9;
	  section.data(10).dtTransOffset = 16;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_tmp
	  section.data(11).logicalSrcIdx = 10;
	  section.data(11).dtTransOffset = 18;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter1_tmp
	  section.data(12).logicalSrcIdx = 11;
	  section.data(12).dtTransOffset = 19;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_tmp_c
	  section.data(13).logicalSrcIdx = 12;
	  section.data(13).dtTransOffset = 20;
	
	  ;% cartpole_MOE_2019a_DW.DiscreteFilter_tmp_a
	  section.data(14).logicalSrcIdx = 13;
	  section.data(14).dtTransOffset = 21;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(1) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.HILInitialize_Card
	  section.data(1).logicalSrcIdx = 14;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(2) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.HILReadEncoderTimebase_Task
	  section.data(1).logicalSrcIdx = 15;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(3) = section;
      clear section
      
      section.nData     = 11;
      section.data(11)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.Force_PWORK.LoggedData
	  section.data(1).logicalSrcIdx = 16;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_DW.HILWriteAnalog_PWORK
	  section.data(2).logicalSrcIdx = 17;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_DW.ToFile_PWORK.FilePtr
	  section.data(3).logicalSrcIdx = 18;
	  section.data(3).dtTransOffset = 2;
	
	  ;% cartpole_MOE_2019a_DW.Voltage_PWORK.LoggedData
	  section.data(4).logicalSrcIdx = 19;
	  section.data(4).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_DW.omega_PWORK.LoggedData
	  section.data(5).logicalSrcIdx = 20;
	  section.data(5).dtTransOffset = 4;
	
	  ;% cartpole_MOE_2019a_DW.theta_PWORK.LoggedData
	  section.data(6).logicalSrcIdx = 21;
	  section.data(6).dtTransOffset = 5;
	
	  ;% cartpole_MOE_2019a_DW.theta_deg_PWORK.LoggedData
	  section.data(7).logicalSrcIdx = 22;
	  section.data(7).dtTransOffset = 6;
	
	  ;% cartpole_MOE_2019a_DW.v_PWORK.LoggedData
	  section.data(8).logicalSrcIdx = 23;
	  section.data(8).dtTransOffset = 7;
	
	  ;% cartpole_MOE_2019a_DW.x_PWORK.LoggedData
	  section.data(9).logicalSrcIdx = 24;
	  section.data(9).dtTransOffset = 8;
	
	  ;% cartpole_MOE_2019a_DW.LQRintegrator_PWORK.LoggedData
	  section.data(10).logicalSrcIdx = 25;
	  section.data(10).dtTransOffset = 9;
	
	  ;% cartpole_MOE_2019a_DW.Scope_PWORK.LoggedData
	  section.data(11).logicalSrcIdx = 26;
	  section.data(11).dtTransOffset = 10;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(4) = section;
      clear section
      
      section.nData     = 4;
      section.data(4)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.HILInitialize_ClockModes
	  section.data(1).logicalSrcIdx = 27;
	  section.data(1).dtTransOffset = 0;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_QuadratureModes
	  section.data(2).logicalSrcIdx = 28;
	  section.data(2).dtTransOffset = 1;
	
	  ;% cartpole_MOE_2019a_DW.HILInitialize_InitialEICounts
	  section.data(3).logicalSrcIdx = 29;
	  section.data(3).dtTransOffset = 3;
	
	  ;% cartpole_MOE_2019a_DW.HILReadEncoderTimebase_Buffer
	  section.data(4).logicalSrcIdx = 30;
	  section.data(4).dtTransOffset = 5;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(5) = section;
      clear section
      
      section.nData     = 1;
      section.data(1)  = dumData; %prealloc
      
	  ;% cartpole_MOE_2019a_DW.ToFile_IWORK.Count
	  section.data(1).logicalSrcIdx = 31;
	  section.data(1).dtTransOffset = 0;
	
      nTotData = nTotData + section.nData;
      dworkMap.sections(6) = section;
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


  targMap.checksum0 = 307575652;
  targMap.checksum1 = 1525571376;
  targMap.checksum2 = 2319926166;
  targMap.checksum3 = 818850156;

