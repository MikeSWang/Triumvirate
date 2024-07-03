"""Custom commit parser for `python-semantic-release` in Conventional
Commits format.

"""
import logging
import re
# from typing import TYPE_CHECKING

from git.objects.commit import Commit
from pydantic.dataclasses import dataclass
from semantic_release import (
    CommitParser,
    LevelBump,
    ParseError,
    ParseResult,
    ParsedCommit,
    ParserOptions,
)
from semantic_release.commit_parser.util import breaking_re, parse_paragraphs

# if TYPE_CHECKING:
#     from git.objects.commit import Commit


TAG_TITLES = {
    'API': "API changes",
    'BUILD': "Maintenance",
    'BUG': "Bug fixes",
    'CI': "None",
    'DOC': "Documentation",
    'ENH': "Improvements",
    'FEAT': "Features",
    'FIX': "Bug fixes",
    'MAINT': "Maintenance",
    'MISC': "Miscellaneous",
    'REV': "None",
    'STYLE': "None",
    'TEST': "None",
}

_commit_filter = "|".join(TAG_TITLES)

_logger = logging.getLogger(__name__)


def _logged_parse_error(commit: Commit,
                        error: str) -> ParseError:  # numpydoc ignore=GL08
    _logger.debug(error)
    return ParseError(commit, error=error)


@dataclass
class TRVParserOptions(ParserOptions):  # numpydoc ignore=GL08
    major_tags: tuple[str, ...] = ('API',)
    minor_tags: tuple[str, ...] = ('FEAT',)
    patch_tags: tuple[str, ...] = ('BUILD', 'BUG', 'ENH', 'FIX', 'MAINT')
    allowed_tags: tuple[str, ...] = (
        # major tags
        'API',
        # minor tags
        'FEAT',
        # patch tags
        'BUILD', 'ENH', 'FIX', 'MAINT',
        # other tags
        'CI', 'DOC', 'MISC', 'REL', 'REV', 'STYLE', 'TEST',
    )
    deprecated_tags: tuple[str, ...] = ('BUG',)
    default_level_bump: LevelBump = LevelBump.NO_RELEASE

    def __post_init__(self) -> None:  # numpydoc ignore=GL08
        self.tag_levels = {
            tag: self.default_level_bump for tag in self.allowed_tags
        }
        for tag in self.patch_tags:
            self.tag_levels[tag] = LevelBump.PATCH
        for tag in self.minor_tags:
            self.tag_levels[tag] = LevelBump.MINOR
        for tag in self.major_tags:
            self.tag_levels[tag] = LevelBump.MAJOR


class TRVCommitParser(CommitParser[ParseResult, TRVParserOptions]):
    """Parser for Triumvirate-style commit messages.

    """  # numpydoc ignore=PR01

    # EXT:TODO: Deprecate in lieu of get_default_options().
    parser_options = TRVParserOptions

    def __init__(self, options: TRVParserOptions | None = None) -> None:

        super().__init__(options)

        self.re_parser = re.compile(
            fr"(?P<type>{_commit_filter})?"
            r"(?:\((?P<scope>[^\n]+)\))?"
            r"(?P<break>!)?:\s+"
            r"(?P<subject>[^\n]+):?"
            r"(\n\n(?P<text>.*))?",
            re.DOTALL,
        )

    @staticmethod
    def get_default_options() -> TRVParserOptions:  # numpydoc ignore=GL08
        return TRVParserOptions()

    def parse(self, commit: Commit) -> ParseResult:  # numpydoc ignore=GL08

        message = str(commit.message)

        parsed = self.re_parser.match(message)
        if not parsed:
            return _logged_parse_error(
                commit, f"Unable to parse the given commit message: {message}"
            )

        parsed_type = parsed.group('type')
        parsed_scope = parsed.group('scope')
        parsed_break = parsed.group('break')
        parsed_subject = parsed.group('subject')
        parsed_text = parsed.group('text')

        # Parse tag title and bump level.
        for tag_ in self.options.allowed_tags + self.options.deprecated_tags:
            if parsed_type == tag_:
                tag_title = TAG_TITLES.get(tag_, "None")
                level_bump = self.options.tag_levels.get(
                    tag_, self.options.default_level_bump
                )
                _logger.debug(
                    "Commit %s introduces a level bump of %s.",
                    commit.hexsha, level_bump,
                )
                break
        else:
            # Some commits may not have a tag, e.g. if they belong to a PR
            # that was not squashed (for maintainability). Ignore these.
            tag_title, level_bump = "None", self.options.default_level_bump
            _logger.debug(
                "Commit %s introduces the default level bump of %s.",
                commit.hexsha, level_bump,
            )

        # Parse commit descriptions.
        if not parsed_subject:
            return _logged_parse_error(
                commit, f"Commit has no subject: {message!r}"
            )

        descriptions = parse_paragraphs(parsed_text) if parsed_text else []
        descriptions.insert(0, parsed_subject)  # insert subject at the top

        # Parse breaking changes.
        breaking_descriptions = [
            match.group(1)
            for match in (breaking_re.match(p) for p in descriptions[1:])
            if match
        ]
        if parsed_break or breaking_descriptions:
            tag_ = 'API'
            tag_title = TAG_TITLES.get(tag_, "API changes")
            level_bump = self.options.tag_levels.get(tag_, LevelBump.MAJOR)
            _logger.debug(
                "Commit %s introduces a breaking-change level bump of %s.",
                commit.hexsha, level_bump
            )

        return ParsedCommit(
            type=tag_title,
            scope=parsed_scope,
            descriptions=descriptions,
            breaking_descriptions=breaking_descriptions,
            bump=level_bump,
            commit=commit,
        )
