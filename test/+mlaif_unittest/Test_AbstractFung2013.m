classdef Test_AbstractFung2013 < matlab.unittest.TestCase
    %% TEST_ABSTRACTFUNG2013 supports mlaif.AbstractFung2013, its class hierarchy, and their software ecosystem.
    %  See also: mlaif.ArterialAnatomy, mlaif.ArterialBoundingBox, mlaif.ArterialSegmentation, mlaif.ArterialCenterline,
    %  mlaif.ECIC, mlaif.ICIC, mlaif.MMRFung2013, mlaif.VisionFung2013, mlaif.Fung2013 (factory).
    %  
    %  Created 23-Mar-2022 15:07:12 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/test/+mlaif_unittest.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    properties
        bids
        cacheMat % assign after segmentations, centerlines, registration targets are ready
        fdg
        fdg_avgt
        ho
        ho_avgt
        mediator
        oc
        oc_avgt
        oo
        oo_avgt
        simplemed
 		testObj
        tof
        t1w
        t1w_os
        t1w_rfov
    end
    
    methods (Test)
        function test_ctor(this)
            disp(this)
        end
        function test_afun(this)
            import mlaif.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);

            this.verifyTrue(isfile(this.fdg.fqfn));
            this.verifyTrue(isfile(this.fdg_avgt.fqfn));
            this.verifyTrue(isfile(this.ho.fqfn));
            this.verifyTrue(isfile(this.oc.fqfn));
            this.verifyTrue(isfile(this.oo.fqfn));
        end
        function test_ArterialAnatomy(this)
            aa = mlaif.ArterialAnatomy.createForT1w('bids', this.simplemed);
            disp(aa)
            this.verifyClass(aa, 'mlaif.ECIC');
            aa.anatomy.view(aa.exclusion_init)
            aa.anatomy.view(aa.inclusion_init)
            this.verifyEqual(aa.anatomy.fileprefix, ...
                'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std')
            this.verifyEqual(aa.exclusion_init.fileprefix, ...
                'wmparc_on_sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std_all')
            this.verifyEqual(aa.inclusion_init.fileprefix, ...
                'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std_b10_thrp33_imdilate_bin')
            this.verifyEqual(aa.wmparc.product.fileprefix, ...
                'wmparc_on_sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std')

            aa = mlaif.ArterialAnatomy.createForTof('bids', this.mediator);
            disp(aa)
            this.verifyClass(aa, 'mlaif.ICIC');
            aa.anatomy.view(aa.exclusion_init)
            aa.anatomy.view(aa.inclusion_init)
            this.verifyEqual(aa.anatomy.fileprefix, ...
                'sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std')
            this.verifyEqual(aa.exclusion_init.fileprefix, ...
                'wmparc_on_sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std_6,28,60,170,173,174,175_imdilate_bin_or_obj')
            this.verifyEqual(aa.inclusion_init.fileprefix, ...
                'sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std_b03_thrp33_imdilate_bin')
            this.verifyEqual(aa.wmparc.product.fileprefix, ...
                'wmparc_on_sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std')
        end
        function test_ArterialBoundingBox(this)
            coord2 = size(this.t1w);
            coord1 = floor(coord2/2);
            abb = mlaif.ArterialBoundingBox(this.t1w, coord1, coord2);
            disp(abb)

            e = abb.extract();
            e = e*4;
            i = abb.insert(e);
            i.view(e)
        end
        function test_ArterialBoundingBox1(this)
            coord1 = [76 154 14];
            coord2 = [62 141 117];
            abb = mlaif.ArterialBoundingBox(this.t1w_os, coord1, coord2);
            disp(abb)

            e = abb.extract();
            e = e*4;
            i = abb.insert(e);
            i.view(e)

            disp(this.t1w_os.imagingFormat.hdr.hist)
            disp(e.imagingFormat.hdr.hist)
            disp(i.imagingFormat.hdr.hist)
            fr = mlfourd.ImagingFormatContext2(fullfile(this.bids.anatPath, 'fslroi_62_141_14_14_13_103_las.nii.gz'));
            disp(fr.hdr.hist)
        end
        function test_ArterialBoundingBox_3D(this)
        end
        function test_ArterialBoundingBox_4D(this)
        end
        function test_ArterialBoundingBox_enlarge_box_base(this)
            abb = mlaif.ArterialBoundingBox(this.t1w, [76 154 14], [62 141 117]);
            e = abb.extract();
            e = e*2;
            ic = abb.insert(e);
            abb = abb.enlarge_box_base();
            e1 = abb.extract();
            e1 = e1*2;
            ic.view(e1);
        end
        function test_ArterialCenterline_T1w(this)
            fung2013 = mlaif.Fung2013.createForT1w('bids', this.simplemed, ...
                'coord1', [76 154 14], ...
                'coord2', [62 141 117]);            
            as = mlaif.ArterialSegmentation('fung2013', fung2013);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, 'segmentation', as.product);
            disp(ac)
            ac = ac.build_centerline();
            disp(ac)

            ac.segmentation.view(ac.product .* 8)
            abb = mlaif.ArterialBoundingBox(ac.anatomy, ac.coord1, ac.coord2);
            anat = abb.insert(ac.product .* dipmax(ac.anatomy));
            anat = anat.threshp(95);
            ac.anatomy.view(anat)

            plot3(ac)
            pcshow(ac)
        end
        function test_ArterialCenterline_tof(this)
            %deleteExisting(fullfile(this.bids.anatPath, 'wmparc_on_sub-*tof*'))
            fung2013 = mlaif.Fung2013.createForTof('bids', this.bids, ...
                'coord1', [228 392 37], ...
                'coord2', [358 515 90]);            
            as = mlaif.ArterialSegmentation('fung2013', fung2013);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, 'segmentation', as.product, 'use_cache', false);
            disp(ac)
            tic
            ac = ac.build_centerline();
            fprintf('ArterialCenterline.build_centerline:')
            toc
            disp(ac)

            %ac.segmentation.view(ac.product .* 8)
            abb = mlaif.ArterialBoundingBox(ac.anatomy, ac.coord1, ac.coord2);
            anat = abb.insert(ac.product .* dipmax(ac.anatomy));
            anat = anat.threshp(95);
            %ac.anatomy.view(anat)

            plot3(ac)
            pcshow(ac)
        end
        function test_ArterialCenterline_tof1(this)
            %deleteExisting(fullfile(this.bids.anatPath, 'wmparc_on_sub-*tof*'))
            fung2013 = mlaif.Fung2013.createForTofOnT1w('bids', this.bids, ...
                'coord1', [63 144 113], ...
                'coord2', [103 182 152]);
            as = mlaif.ArterialSegmentation('fung2013', fung2013);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, 'segmentation', as.product, 'use_cache', false);
            disp(ac)
            tic
            ac = ac.build_centerline();
            fprintf('ArterialCenterline.build_centerline:')
            toc
            disp(ac)

            ac.segmentation.view(ac.product .* 8)
            abb = mlaif.ArterialBoundingBox(ac.anatomy, ac.coord1, ac.coord2);
            anat = abb.insert(ac.product .* dipmax(ac.anatomy));
            anat = anat.threshp(95);
            ac.anatomy.view(anat)

            plot3(ac)
            pcshow(ac)
        end
        function test_ArterialInputFunction_T1w(this)
            fung2013 = mlaif.Fung2013.createForT1w('bids', this.simplemed, ...
                'coord1', [76 154 14], ...
                'coord2', [62 141 117]);
            as = mlaif.ArterialSegmentation('fung2013', fung2013, ...
                'use_cache', true);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'use_cache', true);
            ac = ac.build_centerline();

            aif = mlaif.ArterialInputFunction('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'ic_centerline', ac.product, ...
                'pc_centerline', ac.pointCloud_for_centerline(), ...
                'pet_dyn', this.ho, ...
                'use_cache', true);
            disp(aif)
            aif = aif.build_input_function();
            disp(aif)

            aif.pet_earlyt_on_anatomy1.view(aif.ic_centerline1);
            save(aif.pet_earlyt_on_anatomy1);
            save(aif.ic_centerline1);
            save(aif.product);
            plot(aif)
        end
        function test_ArterialInputFunction_tof(this)
            fung2013 = mlaif.Fung2013.createForTof('bids', this.bids, ...
                'coord1', [228 392 37], ...
                'coord2', [358 515 90]); 
            as = mlaif.ArterialSegmentation('fung2013', fung2013, ...
                'use_cache', true);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'use_cache', true);
            ac = ac.build_centerline();

            aif = mlaif.ArterialInputFunction('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'ic_centerline', ac.product, ...
                'pc_centerline', ac.pointCloud_for_centerline(), ...
                'pet_dyn', this.ho, ...
                'imdilate_scale_mm', 0.2604, ...
                'use_cache', false);
            disp(aif)
            aif = aif.build_input_function();
            disp(aif)

            aif.pet_earlyt_on_anatomy1.view(aif.ic_centerline1);
            save(aif.pet_earlyt_on_anatomy1);
            save(aif.ic_centerline1);
            save(aif.product);
            figure; plot(aif.product);
        end
        function test_ArterialInputFunction_tof1(this)
            fung2013 = mlaif.Fung2013.createForTofOnT1w('bids', this.bids, ...
                'coord1', [63 144 113], ...
                'coord2', [103 182 152]);
            as = mlaif.ArterialSegmentation('fung2013', fung2013, ...
                'use_cache', false);
            as.build_segmentation();

            ac = mlaif.ArterialCenterline('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'use_cache', false);
            ac = ac.build_centerline();

            aif = mlaif.ArterialInputFunction('fung2013', fung2013, ...
                'segmentation', as.product, ...
                'ic_centerline', ac.product, ...
                'pc_centerline', ac.pointCloud_for_centerline(), ...
                'pet_dyn', this.ho, ...
                'use_cache', false);
            disp(aif)
            aif = aif.build_input_function();
            disp(aif)

            aif.pet_earlyt_on_anatomy1.view(aif.ic_centerline1);
            save(aif.pet_earlyt_on_anatomy1);
            save(aif.ic_centerline1);
            save(aif.product);
            figure; plot(aif.product);
        end
        function test_ArterialSegmentation_T1w(this)
            fung2013 = mlaif.Fung2013.createForT1w('bids', this.simplemed, ...
                'coord1', [76 154 14], ...
                'coord2', [62 141 117]);
            as = mlaif.ArterialSegmentation('fung2013', fung2013, 'use_cache', false);
            disp(as)

            [e,bb] = as.build_extracted();
            disp(bb)
            disp(as.anatomy.imagingFormat)
            disp(as.anatomy.imagingFormat.hdr.hist)
            disp(e.imagingFormat)
            disp(e.imagingFormat.hdr.hist)
%            as.anatomy.view(e);

            [em,bb] = as.build_extracted_mask();
            disp(bb)
%            as.anatomy.view(as.exclusion_init)
%            as.anatomy.view(as.inclusion_init)
%            as.anatomy.view(em);

            as.build_segmentation();
            as.product.view(as.anatomy);
            patch(as);
            pcshow(as);
        end
        function test_ArterialSegmentation_tof(this)
            fung2013 = mlaif.Fung2013.createForTof('bids', this.bids, ...
                'coord1', [228 392 37], ...
                'coord2', [358 515 90]);
            as = mlaif.ArterialSegmentation('fung2013', fung2013, 'use_cache', false);
            disp(as)

            [e,bb] = as.build_extracted();
            disp(bb)
%            as.anatomy.view(e);

            [em,bb] = as.build_extracted_mask();
            disp(bb)
%            as.anatomy.view(as.exclusion_init)
%            as.anatomy.view(as.inclusion_init)
%            as.anatomy.view(em);

            as.build_segmentation();
            as.product.view(as.anatomy);
            patch(as);
            pcshow(as);
        end
        function test_bids(this)
            disp(this.bids)
            disp(this.mediator)
        end
        function test_bsplineestim(~)
            % Illustrates B-spline curve estimation without knowing parameter values.
            
            % Copyright 2010 Levente Hunyadi
            
            % spline order
            k = 4;
            % knot sequence
            t = [0 0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1];
            % control points (unknown)
            D_0 = [ 0.1993 0.4965 0.6671 0.7085 0.6809 0.2938 0.1071 0.3929 0.5933 0.8099 0.8998 0.8906 ...
                  ; 0.8377 0.8436 0.7617 0.6126 0.212 0.1067 0.3962 0.5249 0.5015 0.3991 0.6477 0.8553 ];
            % points on B-spline curve
            M_0 = bspline_deboor(k,t,D_0,sort(rand(1,500)));
            M = M_0 + 0.01 * randn(size(M_0));  % contaminate with noise
            
            D = bspline_estimate(k,t,M);
            C = bspline_deboor(k,t,D);
            
            % plot control points and spline
            figure;
            hold all;
            plot(D_0(1,:), D_0(2,:), 'g');
            plot(M_0(1,:), M_0(2,:), 'b');
            plot(M(1,:), M(2,:), 'kx');
            plot(D(1,:), D(2,:), 'r');
            plot(C(1,:), C(2,:), 'c');
            legend('true control points', 'original curve', 'noisy data', 'estimated control points', 'estimated curve', ...
                'Location', 'Best');
            hold off;
        end
        function test_call_t1w(this)
            fung2013 = mlaif.Fung2013.createForT1w('bids', this.simplemed, ...
                'coord1', [76 154 14], ...
                'coord2', [62 141 117], ...
                'timesMid', cumsum(this.simplemed.taus()), ...
                'needs_reregistration', false, ...
                'verbose', 2);
            ic = mlfourd.ImagingContext2(fullfile(this.simplemed.scanPath, 'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet.nii.gz'));
            fung2013.call('pet_dyn', ic, 'use_cache', true);
        end
        function test_call_tof(this)
            fung2013 = mlaif.Fung2013.createForTof('bids', this.bids, ...
                'coord1', [228 392 37], ...
                'coord2', [358 473 90]); % try y = 515 | 473
            fung2013.call('pet_dyn', this.ho, 'use_cache', true);
        end
        function test_createFromCoords(this)
            coords = struct( ...
                't1w', struct( ...
                       'left',  struct( ...
                                'coord1', [139 153 14], ...
                                'coord2', [157 141 117]), ...
                       'right', struct( ...
                                'coord1', [76 154 14], ...
                                'coord2', [62 141 117])), ...
                'tof', struct( ...
                       'left',  struct( ...
                                'coord1', [500 386 37], ...
                                'coord2', [388 520 90]), ...
                       'right', struct( ...
                                'coord1', [228 392 37], ...
                                'coord2', [358 473 90])));
            mlaif.Fung2013.createFromCoords('bids', this.bids, 'coords', coords);
        end
        function test_createFromCoords2(this)
            bids_{1} = mlvg.Ccir1211Bids( ...
                'destinationpath', ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108007'));
            bids_{2} = mlvg.Ccir1211Bids( ...
                'destinationpath', ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108287'));
            bids_{3} = mlvg.Ccir1211Bids( ...
                'destinationpath', ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108300'));

            % for manually selecting coords
%             for i = 1:3
%                 disp(bids_{i}.t1w_ic)
%                 disp(bids_{i}.tof_ic)
%             end
%             return

            coords_{1} = struct( ...
                't1w', struct( ...
                       'left',  struct( ...
                                'coord1', [146 153 47], ...
                                'coord2', [140 163 107]), ...
                       'right', struct( ...
                                'coord1', [93 159 1], ...
                                'coord2', [69 169 101])), ...
                'tof', struct( ...
                       'left',  struct( ...
                                'coord1', [482 377 35], ...
                                'coord2', [385 510 65]), ...
                       'right', struct( ...
                                'coord1', [241 393 25], ...r
                                'coord2', [337 517 59])));
            coords_{2} = struct( ...
                't1w', struct( ...
                       'left',  struct( ...
                                'coord1', [147 123 1], ...
                                'coord2', [128 154 120]), ...
                       'right', struct( ...
                                'coord1', [57 126 1], ...
                                'coord2', [73 145 116])), ...
                'tof', struct( ...
                       'left',  struct( ...
                                'coord1', [479 381 29], ...
                                'coord2', [376 498 56]), ...
                       'right', struct( ...
                                'coord1', [227 360 30], ...
                                'coord2', [327 489 56])));
            coords_{3} = struct( ...
                't1w', struct( ...
                       'left',  struct( ...
                                'coord1', [125 182 73], ...
                                'coord2', [147 156 116]), ...
                       'right', struct( ...
                                'coord1', [83 184 68], ...
                                'coord2', [57 151 113])), ...
                'tof', struct( ...
                       'left',  struct( ...
                                'coord1', [485 401 18], ...
                                'coord2', [395 527 53]), ...
                       'right', struct( ...
                                'coord1', [217 401 19], ...
                                'coord2', [307 529 54])));
            for i = 3:-1:1
                mlaif.Fung2013.createFromCoords('bids', bids_{i}, 'coords', coords_{i});
            end
        end
        function test_Fung2013(this)
            fung2013 = mlaif.Fung2013.createForT1w('bids', this.bids, 'verbose', 0);
            disp(fung2013)
            fung2013 = mlaif.Fung2013.createForTof('bids', this.bids, 'verbose', 2);
            disp(fung2013)
        end
        function test_runSimpleFung2013_mpr(this)
            cd('/data/brier/CerebralMetabolismMS/jjlee/1179204_v1')
            this.simplemed = mlpipeline.SimpleMediator(...
                fullfile(pwd, '1179204_v1_FDG.nii.gz'));
            runSimpleFung2013( ...
                '1179204_v1_FDG.nii.gz', ...
                mpr='1179204_v1_mpr.nii.gz', ...
                mra='1179204_v1_mra.nii.gz', ...
                mpr_coords={[0 0 0]+1, [0 0 0]+1}, ...
                mra_coords={}, ...
                taus=[5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)], ... 
                needs_reregistration=true, ...
                verbose=2, ...
                disp_only=false, ...
                use_cache=false)
        end
        function test_runSimpleFung2013_carotid(this)
            cd('/data/brier/CerebralMetabolismMS/jjlee/1179204_v1_gal')
            this.simplemed = mlpipeline.SimpleMediator(...
                fullfile(pwd, '1179204_v1_FDG.nii.gz'));
            mra_coords = { ...
                [304 231 5; 266 266 96]+1, ...
                [121, 233, 6; 155, 272, 93]+1, ...
                [299, 244, 96; 233, 304, 100]+1, ...
                [123, 249, 93; 188, 304, 97]+1, ...
                [254, 281, 100; 228, 329, 117]+1, ...
                [170, 284, 97; 194, 329, 116]+1};
            runSimpleFung2013( ...
                '1179204_v1_FDG.nii.gz', ...
                mpr='1179204_v1_mpr.nii.gz', ...
                mra='1179204_v1_mra.nii.gz', ...
                mpr_coords={}, ...
                mra_coords=mra_coords, ...
                taus=[5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)], ... 
                needs_reregistration=true, ...
                verbose=1, ...
                disp_only=false, ...
                use_cache=false) 

            % N.B.:
            % [207, 252, 97; 227, 277, 123]+1, basilar
            % [298, 214, 64; 215, 243, 81]+1, R proatlantial
            % [123, 224, 61; 215, 243, 81]+1, L proatlantial
            % [253, 237 0; 276, 215, 35]+1, R vert 
            % [189, 242, 0; 167, 221, 37]+1, L vert 
            % [304 231 5; 266 266 96]+1, R cerv carotid
            % [121, 233, 6; 155, 272, 93]+1, L cerv carotid
            % [299, 244, 96; 233, 304, 100]+1, R petrous
            % [123, 249, 93; 188, 304, 97]+1, L petrous 
            % [254, 281, 100; 228, 329, 117]+1, R supra/clinoid
            % [170, 284, 97; 194, 329, 116]+1, L supra/clinoid
        end
        function test_runSimpleFung2013_vertebral(this)
            cd('/data/brier/CerebralMetabolismMS/jjlee/1179204_v1')
            this.simplemed = mlpipeline.SimpleMediator(...
                fullfile(pwd, '1179204_v1_FDG.nii.gz'));
            mra_coords = { ...
                [298, 214, 64; 215, 243, 81]+1, ...
                [123, 224, 61; 215, 243, 81]+1, ...
                [253, 237 0; 276, 215, 35]+1, ...
                [189, 242, 0; 167, 221, 37]+1};
            runSimpleFung2013( ...
                '1179204_v1_FDG.nii.gz', ...
                mpr='1179204_v1_mpr.nii.gz', ...
                mra='1179204_v1_mra.nii.gz', ...
                mpr_coords={}, ...
                mra_coords=mra_coords, ...
                taus=[5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)], ... 
                needs_reregistration=true, ...
                verbose=1, ...
                disp_only=false, ...
                use_cache=false) 

            % N.B.:
            % [207, 252, 97; 227, 277, 123]+1, basilar
            % [298, 214, 64; 215, 243, 81]+1, R proatlantial
            % [123, 224, 61; 215, 243, 81]+1, L proatlantial
            % [253, 237 0; 276, 215, 35]+1, R vert 
            % [189, 242, 0; 167, 221, 37]+1, L vert 
            % [304 231 5; 266 266 96]+1, R cerv carotid
            % [121, 233, 6; 155, 272, 93]+1, L cerv carotid
            % [299, 244, 96; 233, 304, 100]+1, R petrous
            % [123, 249, 93; 188, 304, 97]+1, L petrous 
            % [254, 281, 100; 228, 329, 117]+1, R supra/clinoid
            % [170, 284, 97; 194, 329, 116]+1, L supra/clinoid
        end

        function test_Wmparc(this)
            w = mlsurfer.Wmparc.createFromBids(this.mediator);
            disp(w)
            w.product.view()

            w = mlsurfer.Wmparc.createCoregisteredFromBids(this.mediator, this.tof);
            disp(w)
            w.product.view()
            deleteExisting(w.wmparc.fqfn)
            deleteExisting(w.T1.fqfn)
        end
        function test_ECIC(this)
            aa = mlaif.ECIC('anatomy', this.t1w, 'bids', this.simplemed);
            disp(aa)
            this.verifyEqual(aa.anatomy.fqfn, ...
                fullfile(this.simplemed.anatPath, strcat(this.t1w.fileprefix, '', '.nii.gz')))
            this.verifyTrue(isfile(aa.anatomy.fqfn))
            aa.anatomy.view(this.t1w) % visually inspect to be superimposable

            try
                deriv_fdg = this.bids.prepare_derivatives(this.fdg_avgt);
            catch
            end
            pet = aa.pet_static_on_anatomy(deriv_fdg);
            pet.view(aa.anatomy)
            ana = aa.anatomy_on_pet(deriv_fdg);
            ana.view(this.fdg_avgt)
        end
        function test_ICIC(this)
            aa = mlaif.ICIC('anatomy', this.tof, 'bids', this.bids);
            disp(aa)
            this.verifyEqual(aa.anatomy.fqfn, ...
                fullfile(this.bids.anatPath, strcat(this.tof.fileprefix, '', '.nii.gz')))
            this.verifyTrue(isfile(aa.anatomy.fqfn))
            aa.anatomy.view(this.tof) % visually inspect to be superimposable

            try
                deriv_fdg = this.bids.prepare_derivatives(this.fdg_avgt);
            catch
            end
            pet = aa.pet_static_on_anatomy(deriv_fdg);
            pet.view(aa.anatomy)
            ana = aa.anatomy_on_pet(deriv_fdg);
            ana.view(this.fdg_avgt)
        end
    end
    
    methods (TestClassSetup)
        function setupAbstractFung2013(this)
 			import mlvg.*;

            % specific to CCIR project 1211, legacy *Bids objects
            this.bids = mlvg.Ccir1211Bids( ...
                'destinationpath', ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108293'));

            % specific to CCIR project 1211, *Mediator design pattern
            this.mediator = mlvg.Ccir1211Mediator( ...
                fullfile(getenv('SINGULARITY_HOME'), 'CCIR_01211', 'derivatives', 'sub-108293', 'ses-20210421', 'pet', ...
                'sub-108293_ses-20210421162709_trc-fdg_proc-static_pet.nii.gz'));

            % working for projects w/o BIDS adherence, sharing base classes with *Mediator
            this.simplemed = mlpipeline.SimpleMediator(...
                fullfile(getenv('SINGULARITY_HOME'), 'SimpleProject', 'sub-simple', ...
                'sub-108293_ses-20210421162709_trc-fdg_proc-static_pet.nii.gz'));
            simplereg = mlpipeline.SimpleRegistry.instance();
            simplereg.tracer = 'FDG';
            simplereg.tausMap = containers.Map;
            simplereg.tausMap('FDG') = [5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)];

            this.fdg =      mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet.nii.gz'));
            this.fdg_avgt = mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421162709_trc-fdg_proc-static_pet.nii.gz'));
            this.ho =       mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421152358_trc-ho_proc-dyn_pet.nii.gz'));
            this.ho_avgt =  mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421152358_trc-ho_proc-static_pet.nii.gz'));
            this.oc =       mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421144815_trc-oc_proc-dyn_pet.nii.gz'));
            this.oc_avgt =  mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421144815_trc-oc_proc-static_pet.nii.gz'));
            this.oo =       mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421154248_trc-oo_proc-dyn_pet.nii.gz'));
            this.oo_avgt =  mlfourd.ImagingContext2(fullfile(this.bids.sourcePetPath, 'sub-108293_ses-20210421154248_trc-oo_proc-static_pet.nii.gz'));            
            this.t1w =      mlfourd.ImagingContext2(fullfile(this.bids.sourceAnatPath, 'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS.nii.gz'));         
            this.t1w_os =   mlfourd.ImagingContext2(fullfile(this.bids.anatPath, 'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_orient-std.nii.gz'));
            this.t1w_rfov = mlfourd.ImagingContext2(fullfile(this.bids.anatPath, 'sub-108293_ses-20210218081506_T1w_MPR_vNav_4e_RMS_robustfov.nii.gz'));
            this.tof =      mlfourd.ImagingContext2(fullfile(this.bids.sourceAnatPath, 'sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab.nii.gz'));
            this.cacheMat = fullfile(getenv('HOME'), 'Tmp', 'Test_AbstractFung2013_cache.mat');
            this.testObj_ = [];
        end
    end
    
    methods (TestMethodSetup)
        function setupAbstractFung2013Test(this)
            if ishandle(this.testObj_)
                this.testObj = copy(this.testObj_);
            end
            %cd(this.bids.anatPath)
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
