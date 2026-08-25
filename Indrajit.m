%% ===== PART 1: Load geometry & identify faces =====
clc; clear; close all;

model = createpde("structural","static-solid");
g = importGeometry(model,"stess ana.stp");

% ---- Step 1: Check actual imported dimensions ----
Vmin = min(g.Vertices);
Vmax = max(g.Vertices);
dims = Vmax - Vmin;
fprintf("Bounding box size: X=%.5f  Y=%.5f  Z=%.5f\n", dims(1), dims(2), dims(3));

figure;
pdegplot(model,"FaceLabels","on","FaceAlpha",0.3);
view(30,25);
title("Check face numbers AND confirm dimensions above match 15x15x2 mm");
%% ===== PART 2 =====
fixedFace = 1;    % from face plot
loadFace  = 3;
loadDir   = [0; -1; 0];

E  = 200e3;   % adjust units to match geometry scale (see note below)
nu = 0.30;
Sy = 250;

structuralProperties(model,"YoungsModulus",E,"PoissonsRatio",nu);
structuralBC(model,"Face",fixedFace,"Constraint","fixed");

A = 15*2
P = 120/A;
structuralBoundaryLoad(model,"Face",loadFace,"SurfaceTraction",P*loadDir);

% ---- Mesh sized relative to actual geometry, not fixed numbers ----
diagLen = norm(dims);          % diagonal of bounding box, real units
Hmax = diagLen/25;             % ~25 elements across the part
Hmin = diagLen/150;            % finer near small features like the bore

generateMesh(model,"Hmax",Hmax,"Hmin",Hmin);

R = solve(model);

figure;
pdeplot3D(model,"ColorMapData",R.VonMisesStress,"Deformation",R.Displacement);
title("Von-Mises Stress"); colorbar;

maxStress = max(R.VonMisesStress);
fprintf("Max Von-Mises Stress = %.4f\n", maxStress);
fprintf("Safety Factor = %.2f\n", Sy/maxStress);

%% ===== Results =====
figure;
pdeplot3D(model,"ColorMapData",R.VonMisesStress,"Deformation",R.Displacement);
title("Von-Mises Stress (MPa)"); colorbar;

maxStress = max(R.VonMisesStress);
maxDisp   = max(sqrt(sum(R.Displacement.ux.^2 + R.Displacement.uy.^2 + R.Displacement.uz.^2,1)));

fprintf("Max Von-Mises Stress = %.2f MPa\n", maxStress);
fprintf("Max Displacement     = %.4f mm\n", maxDisp);
fprintf("Safety Factor        = %.2f\n", Sy/maxStress);