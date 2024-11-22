function [curveData,xp,spatialMesh,params] = readMatcont(system,curve,dim,matContDir,qsave,fileName)

 %Use from the vegPat245 root directory
 %System for Klausmeier: %KlausemeierPTW 
 %system, curve, fileName, matContDir args to be in ''  %marks

 %This script only reads off periodic orbits

 %dim (dimension of the system, should be 3)

 %qsave == 1 to save the data to matContData folder
 %qsave == 0 to not save

 %fileName: if save == 1, give a name to the folder to save the data in

 %matContDir: full path name for where matCont is installed (should end in \matCont7p5)

curveData = load(fullfile(matContDir,'Systems',system,'diagram',curve));
x = curveData.x;

    lds = curveData.globals.lds;
    profile = x(1:lds.ncoords,:);
    xp = reshape(profile,dim,lds.ncol*lds.ntst+1,[]); %xp is a 3d tensor with different variables, points and paramvals

    %f contains the coarse spatial mesh
    f = curveData.f; 
    
    %Remove the last dim rows of f for the spatial mesh for each PO (this is redundant)
    fSize = size(f);
    f = f(1:fSize(1)-dim,:);

    %Matcont algorithm uses ncol evenly spaced points between each entry in f for
    %periodic orbit mesh
    ncol = lds.ncol;

    %Initialise spatialMesh
    spatialMesh = zeros(ncol*(fSize(1)-dim-1)+1,fSize(2));


    %Inserting ncol evenly spaced points between every entry in f
    for j = 1:fSize(2)
        for i = 1:fSize(1)-dim-1
            evenSpacedVector = evenSpacedInsert(f(i,j),f(i+1,j),ncol-1);
            evenSpacedVector = evenSpacedVector(1:ncol); %Remove final entry of this to avoid repetition in spatialMesh
            spatialMesh(ncol*(i-1)+1:ncol*i,j) = evenSpacedVector;
        end
        spatialMesh(ncol*(fSize(1)-dim-1)+1,j) = f(fSize(1)-dim,j);
    end
    %Multiply spatialMesh by the wave speed and period to get the wavelength
    params = x(lds.pars,:);   

    %Save the file if needed
    if qsave == 1
        filePath = fullfile('matContData',fileName);
        save(filePath,'curveData','xp','spatialMesh','params','-mat')
    end
end


