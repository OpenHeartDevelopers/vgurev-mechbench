% 3xn -> 3xn
function f = xyz_to_fiber(xyz)

top = 5;

outer_long_axis = 20;
outer_short_axis = 10;
wall_thickness = 3;
inner_short_axis =  outer_short_axis - wall_thickness;
inner_long_axis =  outer_long_axis - wall_thickness;


% fiber angles
fiber_angle_epi  = -90;
fiber_angle_endo = 90;

x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);

vv=atan2(-y,-x);
%q=(7+3t)  sin u
%z=(17+3t) cos u
%q/(7+3t)=  sin u
%z/(17+3t)= cos u
% q.^2 ./ (7+3*t).^2 + z.^2./(17+3*t).^2 = 1
% (17+3*t).^2 * (7+3*t).^2 - q.^2 .* (17+3*t).^2 - z.^2*(7+3*t).^2 = 0
I = abs(cos(vv)) > 1e-6;
q(I)  = x(I)./cos(vv(I));
q(~I) = y(~I)./sin(vv(~I));

fn = @(t) (inner_long_axis+wall_thickness*t).^2 .* (inner_short_axis+wall_thickness*t).^2 - q.^2 .* (inner_long_axis+wall_thickness*t).^2 - z.^2 .* (inner_short_axis+wall_thickness*t).^2;
t = fsolve(fn,0.5*ones(size(x))); % 4th degree polynomial, fsolve rather than closed-form is used here
uu = acos(z ./ (inner_long_axis+wall_thickness*t));

uu(uu > 0) = -uu(uu > 0); % fix range -pi,0

assert(all(t<1+1e-6 & t>-1e-6)); % check in range
assert(all(uu >= -pi-1e-6));

long_rr  = inner_long_axis  + t*wall_thickness;
short_rr = inner_short_axis + t*wall_thickness;
  

f=zeros(size(xyz));
for i=1:size(uu,2)
  u = uu(i);
  v = vv(i);
  short_r = short_rr(i);
  long_r = long_rr(i);
  fiber_angle = (fiber_angle_endo + t(i) * (fiber_angle_epi - fiber_angle_endo)) * (pi/180);
  deriv_dir = [sin(fiber_angle); cos(fiber_angle)];
  
  M = [short_r.*cos(u)*cos(v) -short_r.*sin(u).*sin(v);  % these are simply the d/du and d/dv in a matrix
       short_r.*cos(u)*sin(v)  short_r.*sin(u).*cos(v);
      -long_r.*sin(u)          0                   ];
  M(:,1) = M(:,1) / norm(M(:,1)); % normalize directions
  M(:,2) = M(:,2) / norm(M(:,2));
  f(:,i) = M * deriv_dir;
  f(:,i) = f(:,i) / norm(f(:,i));
end
quiver3()




