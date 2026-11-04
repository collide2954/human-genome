import js from '@eslint/js';
import globals from 'globals';

const noComments = {
  meta: {
    type: 'problem',
    docs: { description: 'Disallow comments in code' },
    messages: { found: 'Comments are not allowed.' },
    schema: [],
  },
  create(context) {
    return {
      Program() {
        for (const comment of context.sourceCode.getAllComments()) {
          context.report({ loc: comment.loc, messageId: 'found' });
        }
      },
    };
  },
};

export default [
  { ignores: ['node_modules/', 'test-results/'] },
  js.configs.recommended,
  {
    files: ['**/*.{js,mjs}'],
    languageOptions: {
      ecmaVersion: 'latest',
      sourceType: 'module',
      globals: { ...globals.node },
    },
    linterOptions: { reportUnusedDisableDirectives: 'error' },
    plugins: { local: { rules: { 'no-comments': noComments } } },
    rules: {
      'local/no-comments': 'error',
      eqeqeq: 'error',
      'no-var': 'error',
      'prefer-const': 'error',
      'no-shadow': 'error',
      'object-shorthand': 'error',
      'prefer-template': 'error',
      'no-param-reassign': ['error', { props: false }],
      curly: ['error', 'multi-line'],
    },
  },
  {
    files: ['src/**/*.js', 'tests/e2e/**/*.mjs'],
    languageOptions: { globals: { ...globals.browser } },
  },
];
