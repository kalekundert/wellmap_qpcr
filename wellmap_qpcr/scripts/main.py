#!/usr/bin/env python3

import os, sys
import appcli
import matplotlib.pyplot as plt
from pathlib import Path

class App(appcli.App):
    __config__ = [
            appcli.DocoptConfig(),
    ]
    layout_toml = appcli.param('<toml>', cast=Path)
    output = appcli.param('--output', default=None)

    @classmethod
    def main(cls):
        self = cls.from_params()
        appcli.load(self)

        if not self.output:
            if os.fork() != 0:
                sys.exit()

        df, extras = self.load()
        fig = self.plot(df, extras)

        if self.output:
            out = self.output.replace('$', self.layout_toml.stem)
            plt.savefig(out)
            plt.close()
        else:
            plt.show()

