import yaml

class Config():
    def __init__(self, h):
        for k, v in h.iteritems():
            if (type(v) == dict):
                setattr(self, k, Config(v))
            else:
                setattr(self, k, v)

class Settings:
    config = None

    @staticmethod
    def load(path):
        if Settings.config:
            return False

        with open(path, 'r') as f:
            Settings.config = yaml.load(f)
            for k, v in Settings.config.iteritems():
                setattr(Settings, k, Config(v))

        return True
