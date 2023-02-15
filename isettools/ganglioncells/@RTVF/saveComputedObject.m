function RTVFfileWasUpdated = saveComputedObject(obj, computeLconeCenterComputeStruct, computeMconeCenterComputeStruct)
    
   RTVFfileWasUpdated = false;
   if (writeTheFile(obj, computeLconeCenterComputeStruct, computeMconeCenterComputeStruct))
      fprintf('Saving computed object to %s\n', obj.computedObjDataFileName);
      save(obj.computedObjDataFileName, 'obj');
      RTVFfileWasUpdated = true;
   else
      fprintf('Computed RTVF object was not saved to the disk');
   end
   
end

function writeFile = writeTheFile(obj, computeLconeCenterComputeStruct, computeMconeCenterComputeStruct)

   writeFile = false;

   if (~isfile(obj.computedObjDataFileName))
       % Did not find a previously saved file, so we are saving everything
       writeFile = true;
       return;
   end

   % Load the previously saved RTVF file
   fprintf(2,'Found the following previously saved RTVF file\n');
   fprintf(2,'%s\n', obj.computedObjDataFileName);
   theOldOBJ = load(obj.computedObjDataFileName, 'obj');
   
   if (computeLconeCenterComputeStruct) && (computeMconeCenterComputeStruct)
        % Ask user if he wants to overwrite both the L- and the M-cone compute structs
        updateBothLandMstructs = RTVF.queryUserForYesNoResponse('Overwrite both the L- and the M-cone compute structs?');
        if (updateBothLandMstructs)
            writeFile = true;
            return;
        end


        % User wants to overwrite only one of the two. Find out which.
        updateLstruct = RTVF.queryUserForYesNoResponse('Overwrite the L-cone compute struct?');
        updateMstruct = RTVF.queryUserForYesNoResponse('Overwrite the M-cone compute struct?');
        
        if (~updateLstruct) && (~updateMstruct)
           % User  wants to overwrite neither, so do not write the file
           writeFile = false;
           fprintf(2, 'Will not overwrite anything.\n');
           return;
        end
           
        if (updateLstruct) && (~updateMstruct)
           fprintf(2, 'Overwriting the L-cone compute struct, keeping the old M-cone compute struct.\n');
           % Keep the previous M-cone compute struct    
           obj.MconeRFcomputeStruct =  theOldOBJ.obj.MconeRFcomputeStruct; 
           writeFile = true;
           return;
        end
            
        if (~updateLstruct)&&(updateMstruct)
           fprintf(2, 'Overwriting the M-cone compute struct, keeping the old L-cone compute struct.\n');
           % Keep the previous L-cone compute struct    
           obj.LconeRFcomputeStruct =  theOldOBJ.obj.LconeRFcomputeStruct; 
           writeFile = true;
           return;
        end
   end  % (computeLconeCenterComputeStruct) && (computeMconeCenterComputeStruct)


   if (computeLconeCenterComputeStruct) && (~computeMconeCenterComputeStruct)
        updateLstruct = RTVF.queryUserForYesNoResponse('Overwrite the L-cone compute struct?');
        if (updateLstruct)
            fprintf(2, 'Overwriting the L-cone compute struct, keeping the old M-cone compute struct.\n');
            % Keep the previous M-cone compute struct    
            obj.MconeRFcomputeStruct =  theOldOBJ.obj.MconeRFcomputeStruct; 
            writeFile = true;       
        end
        return; 
   end % (computeLconeCenterComputeStruct) && (~computeMconeCenterComputeStruct)

   if (computeMconeCenterComputeStruct) && (~computeLconeCenterComputeStruct)
        updateMstruct = RTVF.queryUserForYesNoResponse('Overwrite the L-cone compute struct?');
        if (updateMstruct)
            fprintf(2, 'Overwriting the M-cone compute struct, keeping the old L-cone compute struct.\n');
            % Keep the previous L-cone compute struct    
            obj.LconeRFcomputeStruct =  theOldOBJ.obj.LconeRFcomputeStruct; 
            writeFile = true;       
        end
        return; 
   end % (computeMconeCenterComputeStruct) && (~computeLconeCenterComputeStruct)

end