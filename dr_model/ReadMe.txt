DR Model
By David Chassin
May 2010

The demand response model provides a time-series model of how electric loads
responds to various demand response signals in an aggregate population.  Three
kinds of demand response controls are handled in this model

1. Direct load control

This mode of demand response involves directly commanding devices to turn off or
on regardless of the prevailing conditions.  In this mode, the fractional rate
at which loads are actuated (on or off) is given per time step.  This mode of 
control is called "ETA" control.

2. Thermostatic reset control

This mode of demand response involves changing the control regime band on the
loads, such as raising or lowering a thermostat setpoint.  In this mode, the 
control band limits are altered.  This mode of control is called "BAND" control.

3. Duty cycle control

This mode of demand response involves modifying the operational duty cycle of
the loads, such as increasing the off time relative to the on time of operation.
This mode of control is called "PHI" control.

All known models of demand response control are some combination of one or more
of these control models and all are supported by the demand response model.

