function [] = outputstate(grid, state)

% Plot some basic fields

subplot(2,3,1)
plot(state.fluid(1).eta,grid.zwkm,'b',state.fluid(2).eta,grid.zwkm,'r')
xlabel('\eta')
ylabel('z')

subplot(2,3,2)
plot(state.fluid(1).q,grid.zwkm,'b',state.fluid(2).q,grid.zwkm,'r')
xlabel('q')
ylabel('z')

subplot(2,3,3)
plot(state.fluid(1).w,grid.zwkm,'b',state.fluid(2).w,grid.zwkm,'r')
xlabel('w')
ylabel('z')

subplot(2,3,4)
plot(state.fluid(1).m,grid.zpkm,'b',state.fluid(2).m,grid.zpkm,'r')
xlabel('m')
ylabel('z')

subplot(2,3,5)
plot(state.p,grid.zpkm,'k')
xlabel('p')
ylabel('z')

subplot(2,3,6)
plot(state.fluid(2).eta-state.fluid(1).eta,grid.zwkm,'k')
xlabel('\Delta \eta')
ylabel('z')

pause(0.01)