% Overwrite selected tendencies with old ones. For checking convergence.
% tend_old needs to be saved in routine advance

disp('*** Overwriting tendencies ***')
% tend.fluid(1).mtke.diffuse = tend_old.fluid(1).mtke.diffuse;
% tend.fluid(1).mtke.diffent = tend_old.fluid(1).mtke.diffent;
% tend.fluid(1).mtke.shear = tend_old.fluid(1).mtke.shear;
% tend.fluid(1).mtke.bflux = tend_old.fluid(1).mtke.bflux;
% tend.fluid(1).mtke.drag = tend_old.fluid(1).mtke.drag;
% tend.fluid(1).mtke.tot = tend.fluid(1).mtke.transport ...
%                       + tend.fluid(1).mtke.diffuse ...
%                       + tend.fluid(1).mtke.diffent ...
%                       + tend.fluid(1).mtke.shear ...
%                       + tend.fluid(1).mtke.bflux ...
%                       + tend.fluid(1).mtke.drag ...
%                       + tend.fluid(1).mtke.dissn ...
%                       + tend.fluid(1).mtke.relabel;
% 
% tend.fluid(2).mtke.transport = tend_old.fluid(2).mtke.transport;
% tend.fluid(2).mtke.diffuse = tend_old.fluid(2).mtke.diffuse;
% tend.fluid(2).mtke.diffent = tend_old.fluid(2).mtke.diffent;
% tend.fluid(2).mtke.dissn   = tend_old.fluid(2).mtke.dissn;
% tend.fluid(2).mtke.relabel = tend_old.fluid(2).mtke.relabel;
% tend.fluid(2).mtke.shear = tend_old.fluid(2).mtke.shear;
% tend.fluid(2).mtke.bflux = tend_old.fluid(2).mtke.bflux;
% tend.fluid(2).mtke.drag = tend_old.fluid(2).mtke.drag;
% tend.fluid(2).mtke.tot = tend.fluid(2).mtke.transport ...
%                        + tend.fluid(2).mtke.diffuse ...
%                        + tend.fluid(2).mtke.diffent ...
%                        + tend.fluid(2).mtke.shear ...
%                        + tend.fluid(2).mtke.bflux ...
%                        + tend.fluid(2).mtke.drag ...
%                        + tend.fluid(2).mtke.dissn ...
%                        + tend.fluid(2).mtke.relabel;
% 
% tend.fluid(1).mvareta.diffuse = tend_old.fluid(1).mvareta.diffuse;
% tend.fluid(1).mvareta.diffent = tend_old.fluid(1).mvareta.diffent;
% tend.fluid(1).mvareta.dissn   = tend_old.fluid(1).mvareta.dissn;
% tend.fluid(1).mvareta.relabel = tend_old.fluid(1).mvareta.relabel;
% tend.fluid(1).mvareta.tot = tend.fluid(1).mvareta.diffuse ...
%                          + tend.fluid(1).mvareta.diffent ...
%                          + tend.fluid(1).mvareta.dissn ...
%                          + tend.fluid(1).mvareta.relabel;
% 
% tend.fluid(2).mvareta.diffuse = tend_old.fluid(2).mvareta.diffuse;
% tend.fluid(2).mvareta.diffent = tend_old.fluid(2).mvareta.diffent;
% tend.fluid(2).mvareta.dissn   = tend_old.fluid(2).mvareta.dissn;
% tend.fluid(2).mvareta.relabel = tend_old.fluid(2).mvareta.relabel;
% tend.fluid(2).mvareta.tot = tend.fluid(2).mvareta.diffuse ...
%                          + tend.fluid(2).mvareta.diffent ...
%                          + tend.fluid(2).mvareta.dissn ...
%                          + tend.fluid(2).mvareta.relabel;
% 
% tend.fluid(2).mw.diffuse = tend_old.fluid(2).mw.diffuse;
% tend.fluid(2).mw.diffent = tend_old.fluid(2).mw.diffent;
% tend.fluid(2).mw.drag    = tend_old.fluid(2).mw.drag;
% tend.fluid(2).mw.pgterm  = tend_old.fluid(2).mw.pgterm;
tend.fluid(2).mw.relabel = tend_old.fluid(2).mw.relabel;
tend.fluid(2).mw.tot = tend.fluid(2).mw.transport ...
                 + tend.fluid(2).mw.diffuse ...
                 + tend.fluid(2).mw.diffent ...
                 + tend.fluid(2).mw.drag ...
                 + tend.fluid(2).mw.pgterm ...
                 + tend.fluid(2).mw.relabel;
% 
% tend.fluid(1).m    = tend_old.fluid(1).m;
% tend.fluid(1).meta = tend_old.fluid(1).meta;
% tend.fluid(1).mq   = tend_old.fluid(1).mq;
% tend.fluid(1).mw   = tend_old.fluid(1).mw;
% tend.fluid(1).mu   = tend_old.fluid(1).mu;
% tend.fluid(1).mv   = tend_old.fluid(1).mv;
% tend.fluid(1).mtke = tend_old.fluid(1).mtke;
% tend.fluid(2).m    = tend_old.fluid(2).m;
% tend.fluid(2).meta = tend_old.fluid(2).meta;
% tend.fluid(2).mq   = tend_old.fluid(2).mq;
% tend.fluid(2).mw   = tend_old.fluid(2).mw;
% tend.fluid(2).mu   = tend_old.fluid(2).mu;
% tend.fluid(2).mv   = tend_old.fluid(2).mv;
% tend.fluid(2).mtke = tend_old.fluid(2).mtke;






