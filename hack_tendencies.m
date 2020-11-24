
% *** Hack tendencies to test budgets ***
if dt <= 2.0
    
    tend.fluid(1).m.relabel = 0*tend.fluid(1).m.relabel;
    tend.fluid(1).m.tot = tend.fluid(1).m.transport ...
                        + tend.fluid(1).m.relabel;
                    
    tend.fluid(1).mw.transport = 0*tend.fluid(1).mw.transport;
    tend.fluid(1).mw.diffuse   = 0*tend.fluid(1).mw.diffuse;
    tend.fluid(1).mw.drag      = 0*tend.fluid(1).mw.drag;
    tend.fluid(1).mw.pgterm    = 1*tend.fluid(1).mw.pgterm;
    tend.fluid(1).mw.relabel   = 1*tend.fluid(1).mw.relabel;
    tend.fluid(1).mw.tot = tend.fluid(1).mw.transport ...
                     + tend.fluid(1).mw.diffuse ...
                     + tend.fluid(1).mw.drag ...
                     + tend.fluid(1).mw.pgterm ...
                     + tend.fluid(1).mw.relabel;
    disp('*** tendencies hacked ***')          
end