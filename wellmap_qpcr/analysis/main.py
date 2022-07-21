#!/usr/bin/env python3

import os, sys
import byoc
import matplotlib.pyplot as plt
from pathlib import Path

class App(byoc.App):
    __config__ = [
            byoc.DocoptConfig,
    ]
    layout_toml = byoc.param('<toml>', cast=Path)
    output = byoc.param('--output', default=None)

    def main(self):
        byoc.load(self)

        if not self.output:
            if os.fork() != 0:
                sys.exit()

        # df, extras = self.load()
        # fig = self.plot(df, extras)

        fig = self.plot(plt.subplots)
        assert fig

        if self.output:
            out = self.output.replace('$', self.layout_toml.stem)
            plt.savefig(out)
            plt.close()
        else:
            plt.show()

