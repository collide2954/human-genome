import stylelint from 'stylelint';

const {
  createPlugin,
  utils: { report, ruleMessages },
} = stylelint;

const ruleName = 'local/no-comments';
const messages = ruleMessages(ruleName, { rejected: 'Comments are not allowed.' });

const rule = (enabled) => (root, result) => {
  if (!enabled) return;
  root.walkComments((node) => {
    report({ ruleName, result, node, message: messages.rejected });
  });
};

rule.ruleName = ruleName;
rule.messages = messages;
rule.meta = { url: 'https://github.com/collide2954/human-genome' };

export default {
  extends: ['stylelint-config-standard'],
  plugins: [createPlugin(ruleName, rule)],
  rules: {
    [ruleName]: true,
    'no-descending-specificity': null,
  },
};
