function saveComputedObject(obj, computeLconeCenterComputeStruct, computeMconeCenterComputeStruct)
    if (isempty(computedRTVObjectExportDirectory))
        fprintf('Did not save the computed object');
        return;
    else
       if (writeTheFile(obj, computeLconeCenterComputeStruct, computeMconeCenterComputeStruct))
           fprintf('Saving computed object to %s\n', obj.computedObjDataFileName);
           save(obj.computedObjDataFileName, 'obj');
       else
           fprintf('Did not save the computed object');
       end
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
   fprintf(2,'Found previously saved RTVF file.\n');
   theOldOBJ = load(obj.computedObjDataFileName, 'obj');
   
   if (computeLconeCenterComputeStruct) && (computeMconeCenterComputeStruct)
        % Ask user if he wants to overwrite both the L- and the M-cone compute structs
        s = input('Overwrite both the L- and the M-cone compute structs? [y=YES] ', 's');
        if (strcmpi(s, 'y'))
            printf(2, 'Overwriting both the L-cone and the M-cone compute struct.\n'); 
            writeFile = true;
            return;
        end

        % User wants to overwrite only one of the two. Find out which.
        overwriteL = false;
        overwriteM = false;
        s = input('Overwrite the L-cone compute struct? [y=YES] ', 's');
        if (strcmpi(s, 'y'))
            overwriteL = true;
        end
        s = input('Overwrite the M-cone compute struct? [y=YES] ', 's');
        if (strcmpi(s, 'y'))
            overwriteM = true;
        end
        
        if (~overwriteL) && (~overwriteM)
           % User  wants to overwrite neither, so do not write the file
           fprintf(2, 'Will not overwrite anything.\n');
           return;
        end
           
        if (overwriteL) && (~overwriteM)
           fprintf(2, 'Overwriting the L-cone compute struct, keeping the old M-cone compute struct.\n');
           % Keep the previous M-cone compute struct    
           obj.MconeRFcomputeStruct =  theOldOBJ.obj.MconeRFcomputeStruct; 
           writeFile = true;
           return;
        end
            
        if (overwriteM)&&(~overwriteL)
           fprintf(2, 'Overwriting the M-cone compute struct, keeping the old L-cone compute struct.\n');
           % Keep the previous L-cone compute struct    
           obj.LconeRFcomputeStruct =  theOldOBJ.obj.LconeRFcomputeStruct; 
           writeFile = true;
           return;
        end
   end  % (computeLconeCenterComputeStruct) && (computeMconeCenterComputeStruct)


   if (computeLconeCenterComputeStruct) && (~computeMconeCenterComputeStruct)
        s = input('Overwrite the L-cone compute struct? [y=YES] ', 's');
        if (strcmpi(s, 'y'))
            fprintf(2, 'Overwriting the L-cone compute struct, keeping the old M-cone compute struct.\n');
            % Keep the previous M-cone compute struct    
            obj.MconeRFcomputeStruct =  theOldOBJ.obj.MconeRFcomputeStruct; 
            writeFile = true;       
        end
        return; 
   end % (computeLconeCenterComputeStruct) && (~computeMconeCenterComputeStruct)

   if (computeMconeCenterComputeStruct) && (~computeLconeCenterComputeStruct)
        s = input('Overwrite the M-cone compute struct? [y=YES] ', 's');
        if (strcmpi(s, 'y'))
            fprintf(2, 'Overwriting the M-cone compute struct, keeping the old L-cone compute struct.\n');
            % Keep the previous L-cone compute struct    
            obj.LconeRFcomputeStruct =  theOldOBJ.obj.LconeRFcomputeStruct; 
            writeFile = true;       
        end
        return; 
   end % (computeMconeCenterComputeStruct) && (~computeLconeCenterComputeStruct)

end