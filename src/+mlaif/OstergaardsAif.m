classdef OstergaardsAif < mlaif.Aif
	%% OSTERGAARDSAIF  
	%  Usage:  obj = OstergaardsAif() 
	%                       ^ 
	%% Version $Revision$ was created $Date$ by $Author$  
	%% and checked into svn repository $URL$ 
	%% Developed on Matlab 7.10.0.499 (R2010a) 
	%% $Id$ 

	properties (Static)
        MaxAuc { get { return maxAuc_; }}
		MinAuc { get { return minAuc_; }}
        maxRoughness_ = Double.NegativeInfinity;
		minRoughness_ = Double.PositiveInfinity;
		maxRoughness { get { return maxRoughness_; }}
		minRoughness { get { return minRoughness_; }}
        
        itsEpiData       { get { return epiData_; }}
		itsPerfusionKeys { get { return pkl_; }}
		itsMasterMask    { get { return masterMask_; }}
		itsConcData      { get { return concData_; }}
		itsAucMap        { get { return aucMap_; }}
		itsRoughnessMap  { get { return roughMap_; }}
		itsClusterMap    { get { return clusterMap_; }}
		itsAifConcCurve  { get { return aifConcCurve_; }}
		itsAifMagnCurve
		{
			get
			{
				ConcModelContext cmc = 
					new ConcModelContext(ref pkl_);
				return cmc.makeMagnCurve(itsAifConcCurve);
			}
		}
		
		clusterDistance
		{
			get { return clusterDistance_; }
			set { clusterDistance_ = value; }
		}
		
		clusterLinkage
		{
			get { return clusterLinkage_; }
			set { clusterLinkage_ = value; }
		}
    end 
    
    properties (Constant)
        FRAC_AUC_TO_TOSS = 0.33;
    end
    
    properties (Static, SetAccess = private)
        epiData_;
		pkl_;
		masterMask_;
		concData_;
		aucMap_;
		roughMap_;
		clusterMap_;
		aifConcCurve_;
		
		clusterDistance_ = Distance.CityBlockFunction;
		clusterLinkage_  = Linkage.CompleteFunction;
    end
    
    properties (SetAccess = private)
        maxAuc_ = Double.NegativeInfinity;
		minAuc_ = Double.PositiveInfinity;
    end

	methods 

		function this = OstergaardsAif(pth, fstem, astem) 
            %% OSTERGAARDSAIF (ctor)  
            %  Usage:  obj = OstergaardsAif() 
            %                       ^ 
            if (nargin < 3); astem = [fstem '_osteraid']; end
            this = this@mlaif.Aif();
        end
        
    end % methods
    
    methods (Static)
        

		

        function cluster = createOstergaardsClustering(IFourdData epidat, ref IPerfusionKeyLists pkl)
		
			epiData_    = new FourdData(epidat);
			pkl_      = pkl;
			concData_   = createConcData(epidat, ref pkl);
			
			initializeMasterMask(concData_.itsKeys);
			
			aucMap_     = createAucMap(concData_);
			updateMasterMask(aucMap_, FRAC_AUC_TO_TOSS*maxAuc_);
			
			roughMap_   = createRoughnessMap(concData_);
			updateMasterMask(roughMap_, FRAC_ROUGHNESS_TO_TOSS*maxRoughness_);

			clusterMap_ = createClusteredMap(concData_, masterMask_); %% Should implement normalization here.
			updateMasterMask(clusterMap_, 0.0);
			
			IConcCurve ccAccum =
				new ConcCurve(new DoubleVector(concData_.lengthT), 
				              pkl.currentConcModel, pkl.currentHctModel, pkl.currentFilterModel);
			for (long v = 0; v < concData_.trueVoxels; v++)
				ccAccum.addAssign(
					new ConcCurve(concData_.evalTimes(Coords4d_red.create(v)), 
								  pkl.currentConcModel, pkl.currentHctModel, pkl.currentFilterModel));
			ccAccum = new ConcCurve(DoubleVector.Divide(ccAccum.itsDoubleVector, clusterMap_.trueVoxels), 
									pkl.currentConcModel, pkl.currentHctModel, pkl.currentFilterModel);
			aifConcCurve_ = new AifConcCurve(ccAccum);
			return aifConcCurve_;
		end
		
		
		
		

		function out = void updateMasterMask(IFourdData map, double lower_limit)
		 
			ICoords4d_red cr;
			ICoords4d     c4d;
			float         ftmp;
			
			if (map.itsFourdStateType == Enums4d.StateType.FOURDFP_RED)
			 
				for (long v = 0; v < map.trueVoxels; v++)
				 
					cr = Coords4d_red.create(v);
					ftmp = (float) map.eval(cr);
					if (ftmp > lower_limit)
						map.assign(1f, cr);
					else map.assign(0f, cr);
				end
			
			else
			 
				for (int z = 0; z < map.lengthZ; z++)
					for (int y = 0; y < map.lengthY; y++)
						for (int x = 0; x < map.lengthX; x++)
						 
							c4d = Coords4d.create(x, y, z);
							ftmp = (float) map.eval(c4d);
							if (ftmp > lower_limit)
								map.assign(1f, c4d);
							else map.assign(0f, c4d);
                            end
                        end
                    end
                end
			end
			map.itsFourdStateType = Enums4d.StateType.FOURDBOOL;
			IFourdbool msk = new Fourdbool((IFourdbool) map.itsFourdState);
			masterMask_ = (IFourdbool) masterMask_.and(msk);
			
			// Also update itsConcData_.itsMask
			concData_ = FourdData.initializeState(concData_.itsFourdState, masterMask_);
		end


		function out =  IFourdData aifConcCurve_to_4dfp(AifConcCurve acc, IInterfileKeys ikeys)
		 
			ikeys.lengths = new int[] { 1, 1, 1, (int) acc.length };
			IFourdData dat = FourdData.createAbInitio(ikeys);
			ICoords4d ct;
			for (int t = 0; t < ikeys.lengthT; t++)
			 
				ct = Coords4d.create(0, 0, 0, t);
				dat.assign((float) acc.itsDoubleVector[t], ct);
			end
			return dat;
		end
		
		

		function out =  IFourdData createConcData(IFourdData epidat, ref IPerfusionKeyLists pkl)
		 
			IFourdData cdata = new FourdData(epidat);
			cdata.itsKeys.comments = "Concentrations by log-fractal strategy";
			cdata.itsKeys.imageModality = "[Gd] by log-fractal strategy";
			cdata.itsKeys.lengthT = pkl.numConcSamples;
			cdata.itsFourdLocator = 
				FourdLocator.updateFilename(epidat.itsFourdLocator, 
				epidat.itsFourdLocator.filenameStem + "_concentrations");
			
			ICoords4d_red ct;
			DoubleVector cVec;
			IConcModel cm = new ConcModelContext(ref pkl);
			for (long v = 0; v < cdata.trueVoxels; v++)
			 
				cVec = cm.makeConcCurve(
						AbsMagnCurve.create(
							 epidat.evalTimes(
							 	Coords4d_red.create(v)))).itsDoubleVector;
				for (int t = 0; t < cdata.lengthT; t++)
				 
					ct = Coords4d_red.create(v, t);
					cdata.assign((float) cVec[t], ct);
				end
			end
			return cdata;
		end
		
		

		function out =  IFourdData createAucMap(IFourdData cdata)
		  
			cdata.trueVoxelsUnknown = true;
			Console.WriteLine("cdata.trueVoxels -> " + cdata.trueVoxels);
			// Console.WriteLine("cdata -> " + cdata);
			
			IFourdData aucmap      = FourdData.initializeState(
												cdata.itsFourdState, cdata.itsMask);
			aucmap.itsKeys         = cdata.itsKeys;
			aucmap.itsKeys.lengthT = 1;
			aucmap.itsFourdLocator =
				FourdLocator.updateFilename(
								cdata.itsFourdLocator, 
								cdata.itsFourdLocator.filenameStem + "_auc");
			aucmap.trueVoxelsUnknown   = true;
			Console.WriteLine("itsAucMap.trueVoxels   -> " + aucmap.trueVoxels);
			// Console.WriteLine("itsAucMap   -> " + itsAucMap);
			
			DoubleVector times = 
				new DoubleVector(cdata.lengthT, 0.0, AbsMagnKeys.T_REPETITION);
			DoubleVector concCurv;
			NaturalCubicSpline concSplines;
			ICoords4d_red cr; float auc;
			for (long v = 0; v < cdata.trueVoxels; v++)
			 
				cr = Coords4d_red.create(v);
				StateAssert.isNotNull(cr, "Coords4d_red object cr was null for v -> " + v);
				concCurv = cdata.evalTimes(cr);
				concSplines = new NaturalCubicSpline(times, concCurv);
				concSplines.Integrator = new GaussKronrodIntegrator();
				auc = (float) concSplines.Integrate(
				              	0.0, cdata.lengthT*AbsMagnKeys.T_REPETITION);
				if (auc > maxAuc_) { maxAuc_ = auc; }
				if (auc < minAuc_) { minAuc_ = auc; }
				try
				 
					aucmap.assign(auc, cr);
				
				catch
				 
					Console.WriteLine("auc -> " + auc + " was unassignable at cr -> " + cr.ToString());
				end
			end
			for (long v = 0; v < aucmap.trueVoxels; v++)
			 
				cr = Coords4d_red.create(v);
				if ((float) aucmap.eval(cr) < FRAC_AUC_TO_TOSS*MaxAuc)
				 
					aucmap.assign(0f, cr);
				end
			end
			aucmap.trueVoxelsUnknown = true;
			aucmap.itsFourdStateType = Enums4d.StateType.FOURDFP_BINARY;
			StateAssert.isTrue(aucmap.trueVoxels <= cdata.trueVoxels, 
				"OstergaardsAifFactory.createAucMap.itsAucMap.trueVoxels was larger than expected -> " + 
				aucmap.trueVoxels);
			return aucmap;
		end
		
		
		
		

		function out =  IFourdData createRoughnessMap(IFourdData cdata)
		 
			IFourdState state    = cdata.itsFourdState; 
			IFourdbool  msk      = new Fourdbool(cdata.itsMask);
			IFourdData  roughmap = FourdData.initializeState(state, msk);
			roughmap.itsKeys     = new InterfileKeys(cdata.itsKeys);
			roughmap.itsKeys.lengthT = 1;
			roughmap.itsFourdLocator =
				FourdLocator.updateFilename(cdata.itsFourdLocator, 
									        cdata.itsFourdLocator.filenameStem + "_roughness");
			
			ICoords4d_red cr;
			double tau = AbsMagnKeys.T_REPETITION;
			DoubleVector times = new DoubleVector(cdata.lengthT, 0.0, tau);
			DoubleVector concCurv;
			DoubleVector firstDeriv, secondDeriv;
			NaturalCubicSpline concSplines, conc1Splines, conc112Splines;
			float roughness;
			roughmap.trueVoxelsUnknown = true;
			cdata.trueVoxelsUnknown     = true;
			Console.WriteLine("cdata.trueVoxels     -> " + cdata.trueVoxels);
			Console.WriteLine("itsRoughnessMap.trueVoxels -> " + roughmap.trueVoxels);
			for (long v = 0; v < roughmap.trueVoxels; v++)
			 
				cr          = Coords4d_red.create(v);
				concCurv    = cdata.evalTimes(cr);
				concSplines = new NaturalCubicSpline(times, concCurv);
				firstDeriv  = new DoubleVector(cdata.lengthT);
				for (int s = 0; s < cdata.lengthT; s++)
					firstDeriv[s] = concSplines.Differentiate(times[s]);
				conc1Splines = new NaturalCubicSpline(times, firstDeriv);
				secondDeriv  = new DoubleVector(cdata.lengthT);
				for (int t = 0; t < cdata.lengthT; t++)
					secondDeriv[t] = conc1Splines.Differentiate(times[t]);
				conc112Splines = new NaturalCubicSpline(times, NMathFunctions.Pow(secondDeriv, 2.0));
				conc112Splines.Integrator = new GaussKronrodIntegrator();
				roughness = (float) conc112Splines.Integrate(0.0, cdata.lengthT*tau);
				if (roughness > maxRoughness_) { maxRoughness_ = roughness; }
				if (roughness < minRoughness_) { minRoughness_ = roughness; }
				// Console.WriteLine("roughness -> " + roughness);
				roughmap.assign(roughness, cr);
			end
			for (long v = 0; v < roughmap.trueVoxels; v++)
			 
				cr = Coords4d_red.create(v);
				if ((float) roughmap.eval(cr) > (1.0 - FRAC_ROUGHNESS_TO_TOSS)*maxRoughness)
				 
					roughmap.assign(0f, cr);
				end
			end
			roughmap.trueVoxelsUnknown = true;
			roughmap.itsFourdStateType = Enums4d.StateType.FOURDFP_BINARY;
			StateAssert.isTrue(roughmap.trueVoxels <= cdata.trueVoxels, 
				"OstergaardsAifFactory.createAucMap.itsAucMap.trueVoxels was larger than expected -> " + roughmap.trueVoxels);
			return roughmap;
		end
		
		
		
		
		
		function out =  IFourdData createNormalizedCurves(IFourdData aucmap, IFourdData cdata)
		 
			IFourdData tmp = new FourdData(cdata);
			ICoords4d_red cvt, cv;
			FloatVector fvec;
			float auc;
			for (long v = 0; v < cdata.trueVoxels; v++)
			 
				cv   = Coords4d_red.create(v);
				fvec = cdata.evalTimes(cv);
				auc  = (float) aucmap.eval(cv);
				if (auc > Math.Sqrt(Single.Epsilon))
					for (long t = 0; t < cdata.lengthT; t++)
					 
						cvt = Coords4d_red.create(v, t);
						tmp.assign((fvec[(int) t])/auc, cvt);
					
				else
					for (long t = 0; t < cdata.lengthT; t++)
					 
						cvt = Coords4d_red.create(v, t);
						tmp.assign(0f, cvt);
					end
			end
			return tmp;
		end
		
		
		
		function out =  IFourdData createClusteredMap(IFourdData cdata, IFourdbool msk)
		 
			cdata = FourdData.initializeState(cdata.itsFourdState, msk);
			IInterfileKeys clusterKeys = new InterfileKeys(cdata.itsKeys);
			clusterKeys.lengthT = 1;
			clusterKeys.numberType = Enums4d.NumberType.FLOAT;
			
			long LENGTH_T = Math.Min(cdata.lengthT, pkl_.maxTimeSamples);
			long MAX_VOXELS = Math.Min(40, cdata.trueVoxels);
			long MIN_VOXELS = 0;
			Console.WriteLine("OstergaardsAifFactory.createClusteredmap.cdata.trueVoxels -> " + cdata.trueVoxels);
			Console.WriteLine("OstergaardsAifFactory.createClusteredmap.MIN_VOXELS          -> " + MIN_VOXELS);
			Console.WriteLine("OstergaardsAifFactory.createClusteredmap.MAX_VOXELS          -> " + MAX_VOXELS);
			
			IFourdData clustMap = FourdData.initializeState(
										Fourdfp_red.createAbInitio(clusterKeys, msk));
			clustMap.itsFourdLocator =
				FourdLocator.updateFilename(cdata.itsFourdLocator, cdata.itsFourdLocator.filenameStem + "_clustered");
			Console.WriteLine("********** created clustermap **************");
			
			DoubleVector   concCurv;
			DoubleMatrix   clusterData   = new DoubleMatrix((int) (MAX_VOXELS - MIN_VOXELS), (int) LENGTH_T);
			ArrayList      clusterAccums = new ArrayList(K_GROUPS);
			long[]         counts        = new long[K_GROUPS];
			
			StateAssert.isNotNull(clusterData, "OstergaardsAifFactory.createClusteredMap.clusterData:  was null (improperly constructed)!");
			Console.WriteLine("*********** assigning clusterAccums **************");
			for (int  k = 0; k < K_GROUPS; k++)
				clusterAccums.Add(new DoubleVector(LENGTH_T));
			Console.WriteLine("*********** assigning counts **************");
			for (int  k = 0; k < K_GROUPS; k++)
				counts[k] = 0;
			Console.WriteLine("*********** assigning concCurv and clusterData **************");
			for (long v = MIN_VOXELS; v < MAX_VOXELS; v++)
			 
				concCurv = cdata.evalTimes(Coords4d_red.create(v));
				StateAssert.isNotNull(concCurv, "OstergaardsAifFactory.createClusterMap.concCurve:  was null at v -> " + v);
				if (v%1000 == 0)
					Console.WriteLine("v -> " + v + " concCurve -> " + Contents.doubleVectorToString(concCurv));
				for (int t = 0; t < cdata.lengthT; t++)
					clusterData[(int) (v - MIN_VOXELS), t] = concCurv[t];
			end
			Console.WriteLine("OstergaardsAifFactory.createClusterMap:   calling ClusterAnalysis ctor............");
			ClusterAnalysis analysis = new ClusterAnalysis(
											clusterData, 
											clusterDistance_, clusterLinkage_);
			StateAssert.isNotNull(analysis, "OstergaardsAifFactory.createClusterMap.analysis:  was null; ClusterAnalysis ctor failed............");
			if (MAX_VOXELS - MIN_VOXELS != analysis.N)
				Messages.showFatalMessage(
					"AifFactor.createClusteredMap.analysis.N was expected to be " + (MAX_VOXELS - MIN_VOXELS) +
					", but was analysis.N -> " + analysis.N);
			Console.WriteLine("OstergaardsAifFactory.createClusterMap:   calling analysis.CutTree............");
			ClusterSet clusterSet = analysis.CutTree(K_GROUPS);
			if (MAX_VOXELS - MIN_VOXELS != clusterSet.N)
				Messages.showFatalMessage(
					"AifFactor.createClusteredMap.clusterSet.N was expected to be " + (MAX_VOXELS - MIN_VOXELS) +
					", but was clusterSet.N -> " + clusterSet.N);
			
			int clust;
			for (long v = MIN_VOXELS; v < MAX_VOXELS; v++)
			 
				clust = clusterSet.Clusters[v];
				if (clust >= K_GROUPS)
					Messages.showFatalMessage(
						"OstergaardsAifFactory.createClusteredMap.clusterSet.Clusters[" + v + 
						"] could not be recognized; clusterSet.Clusters -> " + clusterSet.Clusters[v]);
				clusterAccums[clust] = DoubleVector.Add((DoubleVector) clusterAccums[clust], 
				                                        cdata.evalTimes(Coords4d_red.create(v)));
				counts[clust]++;
			end
			DoubleVector firstMoments = new DoubleVector(K_GROUPS);
			for (int cl = 0; cl < K_GROUPS; cl++)
				DoubleVector.Divide((DoubleVector) clusterAccums[cl], counts[cl]);
			for (int cl = 0; cl < K_GROUPS; cl++)			
				firstMoments[cl] = firstMoment((DoubleVector) clusterAccums[cl]);

			Console.WriteLine("OstergaardsAifFactory.createClusterMap:  assigning itsClusterMap...........");
			int k_best = NMathFunctions.MinIndex(firstMoments);
			Console.WriteLine("************* k_best -> " + k_best + " **************");
			clustMap.assignAll(0f);
			for (long v = 0; v < counts[k_best]; v++)
			 
				clustMap.assign(1f, Coords4d_red.create(v));
			end
			return clustMap;
		end
		
		public const int K_GROUPS = 5;
		
		function out =  double firstMoment(DoubleVector conc)
		 
			DoubleVector times = new DoubleVector(conc.Length, 0.0, AbsMagnKeys.T_REPETITION);
			NaturalCubicSpline integrand = new NaturalCubicSpline(
												times, DoubleVector.Multiply(times, conc));
			integrand.Integrator = new GaussKronrodIntegrator();
			return 
				integrand.Integrate(0.0, conc.Length*AbsMagnKeys.T_REPETITION);
		end
		
		

            end % methods Static





    %  Created with newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end 
