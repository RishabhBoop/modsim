within ModSimBib;
model myDiode
  extends 
    Modelica.Electrical.Analog.Interfaces.OnePort;
protected
  Real s;
equation
  v = s * (if s < 0 then 1 else 0);
  i = s * (if s < 0 then 0 else 1);
end myDiode;