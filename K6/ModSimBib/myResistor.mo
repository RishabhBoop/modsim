within ModSimBib;
model myResistor
  parameter Real R=1;
  extends 
Modelica.Electrical.Analog.Interfaces.OnePort;
equation
  v = R*i;
end myResistor;