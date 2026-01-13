model A2  
  parameter Real pi=Modelica.Constants.pi;
  Real u_e, i_L(start=0), i_e, u_a(start=0);
  parameter Real u_dach=2;
  parameter Real f=1;
  parameter Real R=20;
  parameter Real L=9e-3;
  parameter Real C=1000e-6;
equation
  u_e=u_dach*sin(2*pi*f*time);
  -u_e + R*i_e + u_a=0;
  i_e - C*der(u_a) - i_L=0;
  der(i_L)=u_a/L;
  annotation(
    uses(Modelica(version = "3.2.2")));
end A2;