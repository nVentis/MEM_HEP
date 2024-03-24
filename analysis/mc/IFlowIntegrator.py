from analysis.mc.interfaces import ImportanceSamplingIntegrator

class IFlowIntegrator(ImportanceSamplingIntegrator):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
    