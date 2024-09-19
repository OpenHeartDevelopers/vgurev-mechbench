% Generates cubic lagrange hexahedral meshes for the mechanics benchmark problems 2 and 3
% - cmgui .exnode/.exelem
% - simple text output:
%     .X for nodes
%     .T  for elements (zero-based node indices, 8 corner nodes first in xi-based order followed by other nodes in xi-based order)
%     .T2 for elements (zero-based node indices in xi-based order per element)
%     .F for fibers:  angle x y z per line in the same order as .X
% ".in" problem file for my own code, which may be useful to extract the dirichlet and pressure boundary conditions
%
% Input:
% - nab:   number of elements apex-to-base
% - ncirc: number of elements circumferentially
% - nr:    number of elements transmurally
function [Cnodes,Celems,fib]=benchmark_ellipse_cubic(nab,ncirc,nr)

if nargin<3; nr=2;    end
if nargin<2; ncirc=8; end
if nargin<1; nab=8;   end

output_tag = sprintf('ellipse_benchmark_%d-%d-%d',nr,ncirc,nab);

% ellipse geometry
top = 5;

outer_long_axis = 20;
outer_short_axis = 10;
wall_thickness = 3;

% fiber angles
fiber_angle_epi  = -90;
fiber_angle_endo = 90;


% mesh generation: points to base nodes off

points_circ = ncirc * 3;
points_ab = nab * 3 + 1;
points_trans = nr * 3 + 1;

iz = 1;
coords = zeros(5,points_ab,points_circ+1,iz);

for long_r=(outer_long_axis-wall_thickness):(wall_thickness/(points_trans-1)):outer_long_axis
  short_r = outer_short_axis - (outer_long_axis-long_r);
  
  uu = -acos(top/long_r);
  u = -pi:(pi+uu)/(points_ab-1):uu;
  v = -pi:(2*pi/points_circ):pi;
  
  coords(1,:,:,iz) = short_r .* sin(u)' * cos(v);
  coords(2,:,:,iz) = short_r .* sin(u)' * sin(v);
  coords(3,:,:,iz) = long_r .* cos(u)' * ones(1,length(v));
  coords(4,:,:,iz) = repmat(u',1,length(v));
  coords(5,:,:,iz) = repmat(v,length(u),1);
  
  iz = iz+1;
end

% mesh generation: number nodes naively first, then connect later

el_nodes = zeros(64,5);
naive_elements = zeros(nab*ncirc*nr,64);
naive_nodes    = zeros(64*nab*ncirc*nr,5);
naive_epi1endo0 = zeros(64*nab*ncirc*nr,1);

ei = 0;
for el_ab = 1:nab
  for el_circ = 1:ncirc
    for el_r = 1:nr
      
      ix_ab = (el_ab-1)*3   + (1:4);
      ix_c  = (el_circ-1)*3 + (1:4);
      ix_r  = (el_r-1)*3    + (1:4);
      
      for xi3=1:4
        for xi2=1:4
          for xi1=1:4
            el_nodes(16*(xi3-1) + 4*(xi2-1) + xi1 ,:) = coords(:,ix_ab(xi2),ix_c(xi1),ix_r(xi3));
          end
        end
      end
      
      ei=ei+1;
      nn = (ei-1)*64 + (1:64);
      naive_elements(ei,:) = nn;
      naive_nodes(nn,:) = el_nodes;
      naive_epi1endo0(nn) = (el_r-1)/nr + (floor((0:63)/16)/3)/nr;
    end
  end
end

% merge nodes to create mesh
corner_nodes=[1,4,13,16,49,52,61,64];  no64 = [corner_nodes setdiff(1:64,corner_nodes)];
[Celems,Cnodes,nodemap,revmap] = merge_nodes(naive_elements,naive_nodes(:,1:3));
UV = naive_nodes(revmap,4:5);

% determine angles in actual nodes and check they are consistent
epi1endo0 = zeros(size(Cnodes,1),1);
for i=1:size(Cnodes,1)
  f = naive_epi1endo0(nodemap==i);
  assert( max(f) - min(f) < 1e-10);
  epi1endo0(i) = f(1);
end

% renumbers nodes such that corners are lowest, defining a valid mesh for linear hexes as well
lin = unique(Celems(:,corner_nodes));
to_renumber   = lin(lin > length(lin));
to_renumber_to= setdiff(1:length(lin),lin);

for i=1:length(to_renumber)
  I1 = Celems == to_renumber(i);
  I2 = Celems == to_renumber_to(i);
  Celems(I1) = to_renumber_to(i);
  Celems(I2) = to_renumber(i);
  Cnodes([to_renumber(i) to_renumber_to(i)],:) = Cnodes([to_renumber_to(i) to_renumber(i)],:);
  epi1endo0([to_renumber(i) to_renumber_to(i)]) = epi1endo0([to_renumber_to(i) to_renumber(i)]);
  UV([to_renumber(i) to_renumber_to(i)],:) = UV([to_renumber_to(i) to_renumber(i)],:);
end


fib = zeros(size(Cnodes,1),3);
angle = zeros(size(Cnodes,1),1);
for ni=1:size(Cnodes,1)
  u = UV(ni,1); v = UV(ni,2);
  fiber_angle = (fiber_angle_endo + epi1endo0(ni) * (fiber_angle_epi - fiber_angle_endo)) * (pi/180);
  
  long_r  = outer_long_axis - (1-epi1endo0(ni))*wall_thickness;
  short_r = outer_short_axis - (outer_long_axis-long_r);
  
  deriv_dir = [sin(fiber_angle); cos(fiber_angle)];
  
  M = [short_r*cos(u)*cos(v) -short_r*sin(u)*sin(v);  % these are simply the d/du and d/dv in a matrix
    short_r*cos(u)*sin(v)  short_r*sin(u)*cos(v);
    -long_r*sin(u)          0                   ];
  M(:,1) = M(:,1) / norm(M(:,1)); % normalize directions
  M(:,2) = M(:,2) / norm(M(:,2));
  
  fib(ni,:) = M * deriv_dir;
  fib(ni,:) = fib(ni,:) / norm(fib(ni,:));
  
  angle(ni) = fiber_angle*(180/pi);
  if abs(u+pi)<1e-6 % apex - change this if you want another choice for apical fiber direction
    fib(ni,:) = 0;
    angle(ni) = 0;
  end
end

% ---------- MESH OUTPUT

fh  =fopen(sprintf('out/%s.exelem',output_tag),'w');
fhc =fopen(sprintf('out/%s.T', output_tag),'w');
fhc2=fopen(sprintf('out/%s.T2',output_tag),'w');

fprintf(fh,'Group name: solution\nShape. Dimension=3, line*line*line\n #Scale factor sets=0\n #Nodes=64\n #Fields=1\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,['  ' x '. c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.\n  #Nodes=64\n']);
  for i=1:64
    fprintf(fh,['  ' int2str(i) '. #Values=1\n     Value indices: 1\n     Scale factor indices: 0\n']);
  end
end
fprintf(fhc,'%d\t%d\n',size(Celems,1),size(Cnodes,1));
for e=1:size(Celems,1)
  fprintf(fh,'\nElement: %d 0 0\nNodes: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n',e,Celems(e,:));
  fprintf(fhc ,'%d \n',Celems(e,no64)); fprintf(fhc, '\n');
  fprintf(fhc2,'%d \n',Celems(e,:));    fprintf(fhc2,'\n');
end
fclose(fh);
fclose(fhc);
fclose(fhc2);

fh=fopen(sprintf('out/%s.exnode',output_tag),'w');
fhc=fopen(sprintf('out/%s.X',output_tag),'w');

fprintf(fh,'Region: /\n #Fields=2\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,'  %s.  Value index=%d, #Derivatives=0, #Versions=1\n',x,x-'x'+1);
end
fprintf(fh,'2) fibers, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,'   %s.  Value index=%d, #Derivatives=0, #Versions=1\n',x,x-'x'+4);
end

fprintf(fhc,'%d\t3\n',size(Cnodes,1));
for ni=1:size(Cnodes,1)
  fprintf(fh,'\nNode: %d\n %.8f %.8f %.8f\t%.8f %.8f %.8f',ni,Cnodes(ni,:),fib(ni,:));
  fprintf(fhc,'%.8f %.8f %.8f\n',Cnodes(ni,:));
end
fclose(fh);  fclose(fhc);
fprintf('Wrote exnode/exelem and T/T2/X files for a cubic lagrange hexadral mesh to ./out/%s.*\n',output_tag);


fh=fopen(sprintf('out/%s.in',output_tag),'w');
fprintf(fh,'option: -loadinc 0.001\noption: -inflateonly\n');
fprintf(fh,'mesh: %s\n',output_tag);

topI = find(abs(Cnodes(:,3)-top) < 1e-3);

fprintf(fh,'fixed_nodes: %d\n',length(topI)*3);
for l='xyz'
  fprintf(fh,'%s ',l);
  fprintf(fh,'%d ',topI);
  fprintf(fh,'\n');
end

trac = 1:nr:size(Celems,1);
fprintf(fh,'pressure_boundary 1 %d\n',length(trac));
fprintf(fh,'lv xi3=0 ');
fprintf(fh,'%d ',trac);  % lv pressure bc

fprintf(fh,'\nfields: 5\nepi1endo0 ');
fprintf(fh,'%g ',epi1endo0(1:length(lin)));
fprintf(fh,'\nfiber_angle ');
fprintf(fh,'%g ',angle(1:length(lin)));

for i=1:3
  fprintf(fh,'\nfiber_%s ','x'+(i-1));
  fprintf(fh,'%g ',fib(1:length(lin),i));
end

fprintf('\n');
fclose(fh);
fprintf('Wrote problem file to ./out/%s.in\n',output_tag);

% Merges identical nodes
% nn x 3 -> newnn x 3  . elements size irrelevant
function [elements,nodes,nodemap,reversenodemap] = merge_nodes(elements,nodes)
os = size(nodes,1);
nodes(:,[1 2 3]) = nodes(:,[3 2 1]);
[nodes,reversenodemap,nodemap] = unique(roundn(nodes,-10),'rows'); % with tolerence for rounding
nodes(:,[1 2 3]) = nodes(:,[3 2 1]);

elements = nodemap(elements);
fprintf('removed duplicates %d -> %d\n',os,size(nodes,1));

