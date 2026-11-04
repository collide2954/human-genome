import { Rule, definePlugin } from 'html-validate';

class NoComments extends Rule {
  documentation() {
    return { description: 'Comments are not allowed.' };
  }

  setup() {
    this.on('token', (event) => {
      if (String(event.data[0] ?? '').startsWith('<!--')) {
        this.report(null, 'Comments are not allowed.', event.location);
      }
    });
  }
}

export default {
  extends: ['html-validate:recommended'],
  plugins: [definePlugin({ name: 'local', rules: { 'local/no-comments': NoComments } })],
  rules: { 'local/no-comments': 'error' },
};
